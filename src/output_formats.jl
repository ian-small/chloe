
using BioSequences
using GenomicAnnotations
using UUIDs
using Dates

function write_model2SFF(outfile::IO, model::SFF_Model)
    isnothing(model) && return
    model_id = model.gene * "/" * string(model.gene_count)
    for sff in model.features
        f = sff.feature
        id = "$(model.gene)/$(model.gene_count)/$(f.type)/$(f.order)"
        write(outfile, id)
        write(outfile, "\t")
        write(outfile, join([model.strand, string(f.start), string(f.length), string(f.phase)], "\t"))
        write(outfile, "\t")
        write(
            outfile,
            join(
                [
                    @sprintf("%.3g", sff.relative_length),
                    @sprintf("%.3g", sff.stackdepth),
                    @sprintf("%.3g", sff.gmatch),
                    @sprintf("%.3g", sff.feature_prob),
                    @sprintf("%.3g", sff.coding_prob)
                ],
                "\t"
            )
        )
        write(outfile, "\t")
        write(outfile, join(model.warnings, "; "))
        write(outfile, "\n")
    end
    for warning in model.warnings
        @warn("$(model_id) $(warning)")
    end
end

function writeSFF(
    outfile::Union{String,IO},
    id::String, # NCBI id
    genome_length::Int32,
    mean_coverage::Float32,
    models::FwdRev{Vector{SFF_Model}}
)
    function out(outfile::IO)
        write(outfile, id, "\t", string(genome_length), "\t", @sprintf("%.3f", mean_coverage), "\n")
        for model in models.forward
            write_model2SFF(outfile, model)
        end
        for model in models.reverse
            write_model2SFF(outfile, model)
        end
    end
    if outfile isa String
        maybe_gzwrite(outfile::String) do io
            out(io)
        end
    else
        out(outfile)
    end
end

function construct_locus(sff::SFF_Model, feature_type::String, target_length::Integer)
    sfeatures = filter(x -> x.feature.type == feature_type, sff.features)
    isempty(sfeatures) && return nothing
    loci = Vector{AbstractLocus}()
    for sfeature in sfeatures
        f = sfeature.feature
        start = sff.strand == '+' ? f.start : reverse_complement(f.start + f.length-1, target_length)
        span = ClosedSpan(start:start + f.length-1)
        push!(loci, span)
    end
    locus = length(loci) == 1 ? first(loci) : Join(sff.strand == '+' ? loci : reverse(loci))
    if sff.strand == '-'; locus = Complement(locus); end
    locus
end

function chloe2biojulia(chloe::ChloeAnnotation)::GenomicAnnotations.Record
    biojulia = GenomicAnnotations.Record()
    biojulia.name = chloe.target_id
    biojulia.header = "##gff-version 3\n##source-version Chloe 1.0\n##sequence-region\t$(chloe.target_id)\t1\t$(chloe.target_length)\n"
    biojulia.circular = true
    # add

    # define a product dictionary (not sure this is the right place to do this)
    product_dict = Dict(
        "accD"=>"acetyl-coenzyme A carboxylase carboxyl transferase subunit beta",
        "atpA"=>"ATP synthase subunit alpha",
        "atpB"=>"ATP synthase subunit beta",
        "atpE"=>"ATP synthase epsilon chain",
        "atpF"=>"ATP synthase subunit b",
        "atpH"=>"ATP synthase subunit c",
        "atpI"=>"ATP synthase subunit a",
        "ccsA"=>"cytochrome c biogenesis protein",
        "cemA"=>"potassium/proton antiporter",
        "clpP1"=>"chloroplastic ATP-dependent Clp protease proteolytic subunit 1",
        "infA"=>"translation initiation factor IF-1",
        "matK"=>"maturase K",
        "ndhA"=>"NAD(P)H-quinone oxidoreductase subunit 1",
        "ndhB"=>"NAD(P)H-quinone oxidoreductase subunit 2",
        "ndhC"=>"NAD(P)H-quinone oxidoreductase subunit 3",
        "ndhD"=>"NAD(P)H-quinone oxidoreductase chain 4",
        "ndhE"=>"NAD(P)H-quinone oxidoreductase subunit 4L",
        "ndhF"=>"NAD(P)H-quinone oxidoreductase subunit 5",
        "ndhG"=>"NAD(P)H-quinone oxidoreductase subunit 6",
        "ndhH"=>"NAD(P)H-quinone oxidoreductase subunit H",
        "ndhI"=>"NAD(P)H-quinone oxidoreductase subunit I",
        "ndhJ"=>"NAD(P)H-quinone oxidoreductase subunit J",
        "ndhK"=>"NAD(P)H-quinone oxidoreductase subunit K",
        "pafI"=>"photosystem I assembly protein Ycf3",
        "pafII"=>"photosystem I assembly protein Ycf4",
        "pbf1"=>"photosystem biogenesis factor 1",
        "petA"=>"cytochrome f",
        "petB"=>"cytochrome b6",
        "petD"=>"cytochrome b6-f complex subunit 4 ",
        "petG"=>"cytochrome b6-f complex subunit 5",
        "petL"=>"cytochrome b6-f complex subunit 6",
        "petN"=>"cytochrome b6-f complex subunit 8",
        "psaA"=>"photosystem I P700 chlorophyll a apoprotein A1",
        "psaB"=>"photosystem I P700 chlorophyll a apoprotein A2",
        "psaC"=>"photosystem I iron-sulfur center",
        "psaI"=>"photosystem I reaction center subunit VIII",
        "psaJ"=>"photosystem I reaction center subunit IX",
        "psbA"=>"photosystem II protein D1",
        "psbB"=>"photosystem II CP47 reaction center protein",
        "psbC"=>"photosystem II CP43 reaction center protein",
        "psbD"=>"photosystem II D2 protein",
        "psbE"=>"cytochrome b559 subunit alpha",
        "psbF"=>"cytochrome b559 subunit beta",
        "psbH"=>"photosystem II reaction center protein H",
        "psbI"=>"photosystem II reaction center protein I",
        "psbJ"=>"photosystem II reaction center protein J",
        "psbK"=>"photosystem II reaction center protein K",
        "psbL"=>"photosystem II reaction center protein L",
        "psbM"=>"photosystem II reaction center protein M",
        "psbT"=>"photosystem II reaction center protein T",
        "psbZ"=>"photosystem II reaction center protein Z",
        "rbcL"=>"ribulose bisphosphate carboxylase large chain",
        "rpl14"=>"large ribosomal subunit protein uL14c",
        "rpl16"=>"large ribosomal subunit protein uL16c",
        "rpl2"=>"large ribosomal subunit protein uL2cz/uL2cy",
        "rpl20"=>"large ribosomal subunit protein bL20c",
        "rpl22"=>"large ribosomal subunit protein uL22c",
        "rpl23"=>"large ribosomal subunit protein uL23cz/uL23cy",
        "rpl32"=>"large ribosomal subunit protein bL32c",
        "rpl33"=>"large ribosomal subunit protein bL33c",
        "rpl36"=>"large ribosomal subunit protein bL36c",
        "rpoA"=>"DNA-directed RNA polymerase subunit alpha",
        "rpoB"=>"DNA-directed RNA polymerase subunit beta",
        "rpoC1"=>"DNA-directed RNA polymerase subunit beta'",
        "rpoC2"=>"DNA-directed RNA polymerase subunit beta''",
        "rps11"=>"small ribosomal subunit protein uS11c",
        "rps12"=>"small ribosomal subunit protein uS12cz/uS12cy",
        "rps14"=>"small ribosomal subunit protein uS14c",
        "rps15"=>"small ribosomal subunit protein uS15c",
        "rps16"=>"small ribosomal subunit protein uS16c",
        "rps18"=>"small ribosomal subunit protein uS18c",
        "rps19"=>"small ribosomal subunit protein uS19c",
        "rps2"=>"small ribosomal subunit protein uS2c",
        "rps3"=>"small ribosomal subunit protein uS3c",
        "rps4"=>"small ribosomal subunit protein uS4c",
        "rps7"=>"small ribosomal subunit protein uS7cz/uS7cy",
        "rps8"=>"small ribosomal subunit protein uS8c",
        "ycf1"=>"protein TIC 214",
        "ycf2"=>"protein Ycf2",
        "rrn16"=>"16S ribosomal RNA",
        "rrn23"=>"23S ribosomal RNA",
        "rrn4.5"=>"4.5S ribosomal RNA",
        "rrn5"=>"5S ribosomal RNA",
        "trnA-UGC"=>"tRNA-Ala",
        "trnC-GCA"=>"tRNA-Cys",
        "trnD-GUC"=>"tRNA-Asp",
        "trnE-UUC"=>"tRNA-Glu",
        "trnF-GAA"=>"tRNA-Phe",
        "trnfM-CAU"=>"tRNA-fMet",
        "trnG-GCC"=>"tRNA-Gly",
        "trnG-UCC"=>"tRNA-Gly",
        "trnH-GUG"=>"tRNA-His",
        "trnI-CAU"=>"tRNA-Ile",
        "trnI-GAU"=>"tRNA-Ile",
        "trnK-UUU"=>"tRNA-Lys",
        "trnL-CAA"=>"tRNA-Leu",
        "trnL-UAA"=>"tRNA-Leu",
        "trnL-UAG"=>"tRNA-Leu",
        "trnM-CAU"=>"tRNA-Met",
        "trnN-GUU"=>"tRNA-Asn",
        "trnP-UGG"=>"tRNA-Pro",
        "trnQ-UUG"=>"tRNA-Gln",
        "trnR-ACG"=>"tRNA-Arg",
        "trnR-UCU"=>"tRNA-Arg",
        "trnS-GCU"=>"tRNA-Ser",
        "trnS-GGA"=>"tRNA-Ser",
        "trnS-UGA"=>"tRNA-Ser",
        "trnT-GGU"=>"tRNA-Thr",
        "trnT-UGU"=>"tRNA-Thr",
        "trnV-GAC"=>"tRNA-Val",
        "trnV-UAC"=>"tRNA-Val",
        "trnW-CCA"=>"tRNA-Trp",
        "trnY-GUA"=>"tRNA-Tyr"
    )

    sffs = vcat(chloe.annotation.forward, chloe.annotation.reverse)
    locus_index = 1
    for sff in sffs
        merge_adjacent_features!(sff)
        startswith(sff.gene, "rps12") && continue
        ft = featuretype(sff)
        if ft ≠ "repeat_region"; ft = "gene"; end
        # construct gene/repeat_region feature
        span = gene_span(sff)
        if sff.strand == '-'; span = reverse_complement(span, chloe.target_length); end
        locus = ClosedSpan(span)
        if sff.strand == '-'; locus = Complement(locus); end
        gene_id = string(uuid4())
        locus_tag = "LOCUSTAG_" * string(locus_index) * "loc"
        addgene!(biojulia, Symbol(ft), locus; locus_tag = locus_tag, ID = gene_id, gene = sff.gene, Name = sff.gene)
        # optionally construct mRNA feature
        # construct CDS, tRNA or rRNA feature
        for feature_type in ["CDS", "tRNA", "rRNA"]
            if feature_type == "CDS"
                locus = construct_locus(sff, feature_type, chloe.target_length)
                if ~isnothing(locus)
                    id = string(uuid4())
                    addgene!(biojulia, Symbol(feature_type), locus; Parent = gene_id, locus_tag = locus_tag, ID = id, gene = sff.gene,
                        Name = "$(sff.gene).$feature_type", product = get(product_dict, sff.gene, "missing"), transl_table = 11)
                end
            else
                locus = construct_locus(sff, feature_type, chloe.target_length)
                if ~isnothing(locus)
                    id = string(uuid4())
                    addgene!(biojulia, Symbol(feature_type), locus; Parent = gene_id, locus_tag = locus_tag, ID = id, gene = sff.gene,
                        Name = "$(sff.gene).$feature_type", product = get(product_dict, sff.gene, "missing"))
                end
            end
        end
        # construct intron feature(s)
        introns = filter(x -> x.feature.type == "intron", sff.features)
        for (i, intron) in enumerate(introns)
            f = intron.feature
            start = sff.strand == '+' ? f.start : reverse_complement(f.start + f.length-1, chloe.target_length)
            locus = ClosedSpan(start:start + f.length-1)
            if sff.strand == '-'; locus = Complement(locus); end
            id = string(uuid4())
            addgene!(biojulia, :intron, locus; Parent = gene_id, locus_tag = locus_tag, ID = id, gene = sff.gene, Name = "$(sff.gene).intron.$i",
                number = i)
        end
        locus_index += 1
    end
    # join rps12A and rps12B features
    for a in filter(x -> x.gene == "rps12A", sffs), b in filter(x -> x.gene == "rps12B", sffs)
        locus_tag = "LOCUSTAG_" * string(locus_index) * "loc"
        # construct gene/repeat_region feature
        aspan = gene_span(a)
        if a.strand == '-'; aspan = reverse_complement(aspan, chloe.target_length); end
        alocus = ClosedSpan(aspan)
        if a.strand == '-'; alocus = Complement(alocus); end
        bspan = gene_span(b)
        if b.strand == '-'; bspan = reverse_complement(bspan, chloe.target_length); end
        blocus = ClosedSpan(bspan)
        if b.strand == '-'; blocus = Complement(blocus); end
        gene_id = string(uuid4())
        addgene!(biojulia, :gene, Join([alocus, blocus]); locus_tag = locus_tag, ID = gene_id, gene = "rps12", Name = "rps12")
        # construct CDS feature
        alocus = construct_locus(a, "CDS", chloe.target_length)
        bloci = Vector{AbstractLocus}()
        for sfeature in filter(x -> x.feature.type == "CDS", b.features)
            f = sfeature.feature
            start = b.strand == '+' ? f.start : reverse_complement(f.start + f.length-1, chloe.target_length)
            span = ClosedSpan(start:start + f.length-1)
            push!(bloci, b.strand == '+' ? span : Complement(span))
        end
        id = string(uuid4())
        addgene!(biojulia, :CDS, Join([alocus, bloci...]); Parent = gene_id, locus_tag = locus_tag, ID = id, gene = "rps12", Name = "rps12.CDS",
            product = "small ribosomal subunit protein uS12cz/uS12cy", transl_table = 11)
        # construct rps12A internal intron feature(s) (I don't think there are any, but just in case...)
        aintrons = filter(x -> x.feature.type == "intron", a.features)
        intron_count = 1
        for intron in aintrons[1:end-1]
            f = intron.feature
            start = a.strand == '+' ? f.start : reverse_complement(f.start + f.length-1, chloe.target_length)
            locus = ClosedSpan(start:start + f.length-1)
            if a.strand == '-'; locus = Complement(locus); end
            id = string(uuid4())
            addgene!(biojulia, :intron, locus; Parent = gene_id, locus_tag = locus_tag, ID = id, gene = "rps12", Name = "rps12.intron.$intron_count",
                number = intron_count)
            intron_count += 1
        end
        bintrons = filter(x -> x.feature.type == "intron", b.features)
        # construct transpliced intron feature
        if ~isempty(aintrons)
            f = last(aintrons).feature
            start = a.strand == '+' ? f.start : reverse_complement(f.start + f.length-1, chloe.target_length)
            alocus = ClosedSpan(start:start + f.length-1)
            if a.strand == '-'; alocus = Complement(alocus); end
            if ~isempty(bintrons)
                f = first(bintrons).feature
                start = b.strand == '+' ? f.start : reverse_complement(f.start + f.length-1, chloe.target_length)
                blocus = ClosedSpan(start:start + f.length-1)
                if b.strand == '-'; blocus = Complement(blocus); end
                id = string(uuid4())
                addgene!(biojulia, :intron, Join([alocus, blocus]); Parent = gene_id, locus_tag = locus_tag, ID = id, gene = "rps12",
                    Name = "rps12.intron.$intron_count", number = intron_count)
                intron_count += 1
            end
        end
        # construct rps12B internal intron feature(s)
        for intron in bintrons[2:end]
            f = intron.feature
            start = b.strand == '+' ? f.start : reverse_complement(f.start + f.length-1, chloe.target_length)
            locus = ClosedSpan(start:start + f.length-1)
            if b.strand == '-'; locus = Complement(locus); end
            id = string(uuid4())
            addgene!(biojulia, :intron, locus; Parent = gene_id, locus_tag = locus_tag, ID = id, gene = "rps12", Name = "rps12.intron.$intron_count",
                number = intron_count)
            intron_count += 1
        end
    end
    sort!(biojulia.genes)
    biojulia
end

function write_result(config::ChloeConfig, target::FwdRev{CircularSequence}, result::ChloeAnnotation, filestem::String)::Tuple{Union{String,IO},String}
    if ~config.no_transform
        FASTAWriter(open(filestem * ".chloe.fa", "w")) do outfile
            write(outfile, FASTARecord(result.target_id, target.forward[1:length(target.forward)]))
        end
    end
    if config.sff
        out = filestem * ".chloe.sff"
        writeSFF(out, result.target_id, result.target_length, geomean(values(result.coverages)), result.annotation)
    end
    if ~config.no_gff || config.gbk || config.embl
        biojulia = chloe2biojulia(result)
        if ~config.no_gff
            biojulia.sequence = dna""
            out = filestem * ".chloe.gff"
            GFF.printgff(out, biojulia)
        end
        if config.gbk
            out = filestem * ".chloe.gbk"
            biojulia.header = "LOCUS       $(rpad(result.target_id, 10, ' ')) $(lpad(length(biojulia.sequence), 10, ' ')) bp    DNA     circular PLN $(uppercase(Dates.format(now(), "dd-uuu-yyyy")))"
            biojulia.sequence = target.forward[1:length(target.forward)]
            GenBank.printgbk(out, biojulia)
        end
        if config.embl
            out = filestem * ".chloe.embl"
            biojulia.sequence = target.forward[1:length(target.forward)]
            EMBL.printembl(out, biojulia)
        end
    end
    return out, result.target_id
end
