
module CmdLine
export chloe_main

import ArgParse: ArgParseSettings, @add_arg_table!, parse_args
import Logging

import ..Annotator

include("../globals.jl")
include("dist_globals.jl")

function quiet_metafmt(level, _module, group, id, file, line)
    color = Logging.default_logcolor(level)
    prefix = (level == Logging.Warn ? "Warning" : string(level)) * ':'
    return color, prefix, ""
end

function chloe(; reference_dir="default", fasta_files=String[],
    sensitivity=DEFAULT_SENSITIVITY,
    output::Union{Nothing,String}=nothing, gff::Bool=false, nofilter::Bool=false)
    db = Annotator.ReferenceDb(reference_dir=reference_dir)
    config = Annotator.ChloeConfig(; sensitivity=sensitivity, to_gff3=gff, nofilter=nofilter)
    Annotator.annotate_batch(db, fasta_files, config, output)
end

function getargs(args::Vector{String}=ARGS)
    cmd_args = ArgParseSettings(prog="Chloë", autofix_names=true)  # turn "-" into "_" for arg names.

    @add_arg_table! cmd_args begin
        "minhash"
        help = "minhash fasta files of reference genomes"
        action = :command
        "align"
        help = "align 2 chloroplast genomes"
        action = :command
        "annotate"
        help = "annotate fasta files"
        action = :command
        "rotate"
        help = "rotate circular genomes to standard position"
        action = :command
        "--level", "-l"
        arg_type = String
        default = "info"
        help = "log level (info,warn,error,debug)"
    end

    @add_arg_table! cmd_args["minhash"] begin
        "fasta-files"
        arg_type = String
        nargs = '+'
        required = true
        action = :store_arg
        help = "fasta files to process (or directory of files)"
    end

    @add_arg_table! cmd_args["align"] begin
        "query"
        arg_type = String
        required = true
        action = :store_arg
        help = "query sequence (fasta format)"
        "target"
        arg_type = String
        nargs = '+'
        required = true
        action = :store_arg
        help = "target sequence(s): fasta file(s) to process (or directory of files)"
        "--output", "-o"
        arg_type = String
        default = "default"
        help = "output file"
    end

    @add_arg_table! cmd_args["annotate"] begin
        "fasta-files"
        arg_type = String
        nargs = '+'
        required = true
        action = :store_arg
        help = "fasta files to process"
        "--output", "-o"
        arg_type = String
        help = "output filename (or directory if multiple fasta files)"
        "--reference", "-r"
        arg_type = String
        default = "cp"
        dest_name = "reference_dir"
        metavar = "DIRECTORY"
        help = "references and templates to use for annotations: cp for chloroplast, nr for nuclear rDNA [default: cp]"
        "--sensitivity", "-s"
        arg_type = Float64
        default = DEFAULT_SENSITIVITY
        help = "probability threshold for reporting features [default: $(DEFAULT_SENSITIVITY)]"
        "--nofilter"
        action = :store_true
        help = "don't filter output"
        "--gff"
        action = :store_true
        help = "save output in gff3 format instead of sff"
    end

    @add_arg_table! cmd_args["rotate"] begin
        "fasta-files"
        arg_type = String
        nargs = '+'
        required = true
        action = :store_arg
        help = "fasta file(s) to process"
        "--flip-SSC", "-S"
        action = :store_true
        help = "flip orientation of small single-copy region"
        "--flip-LSC", "-L"
        action = :store_true
        help = "flip orientation of large single-copy region"
        "--extend", "-e"
        default = 0
        arg_type = Int
        help = "add n bases from start to the end of sequence to allow mapping to wrap ends [use -1 for maximum extent]"
        "--output", "-o"
        arg_type = String
        help = "output file or directory"
    end

    # args.epilog = """
    #     examples:\n
    #     \ua0\ua0 # chloe.jl -t template.tsv -r reference_dir fasta1 fasta2 ...\n
    #     """
    parse_args(args, cmd_args; as_symbols=true)
end

function chloe_main(args::Vector{String} = ARGS)
    parsed_args = getargs(args)
    level = lowercase(parsed_args[:level])
    Logging.with_logger(Logging.ConsoleLogger(stderr,
        get(LOGLEVELS, level, Logging.Warn), meta_formatter=quiet_metafmt)) do
        cmd = parsed_args[:_COMMAND_]
        a = parsed_args[cmd]
        if cmd == :minhash
            Annotator.minhash_references(; a...)
        elseif cmd == :align
            Annotator.align(; a...)
        elseif cmd == :annotate
            chloe(; a...)
        elseif cmd == :rotate
            Annotator.rotategenomes(; a...)
        end
    end

end

end # module

