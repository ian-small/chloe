module Broker

export broker_main, test_bind

import ZMQ
import ZeroMQ_jll
import ArgParse: ArgParseSettings, @add_arg_table!, parse_args

function test_bind(worker_url::String, client_url::String)
    # this doesn't seem to work
    # a test bind to an already bound
    # ipc:///named-socket seems to stuff
    # everything up.....
    router = ZMQ.Socket(ZMQ.ROUTER)
    dealer = ZMQ.Socket(ZMQ.DEALER)
    try
        ZMQ.bind(router, client_url)
        ZMQ.bind(dealer, worker_url)
        return nothing
    catch e
        return e.msg
    finally
        ZMQ.set_linger(router, 0)
        ZMQ.set_linger(dealer, 0)
        ZMQ.close(router)
        ZMQ.close(dealer)
    end

end
function start_broker(worker_url::String, client_url::String)

    # ctx = Context()
    router = ZMQ.Socket(ZMQ.ROUTER)
    dealer = ZMQ.Socket(ZMQ.DEALER)

    ZMQ.bind(router, client_url)
    ZMQ.bind(dealer, worker_url)

    # missing: ZMQ.proxy(router, dealer)
    try
        try
            @info "worker=$worker_url client=$client_url"
            rc = ccall((:zmq_proxy, ZeroMQ_jll.libzmq), Cint,  (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}), router, dealer, C_NULL)
            @info "done zmq_proxy $(rc)"
        catch
            @error "zmq_proxy exception!"
        end
    finally
        # control never comes here... clean up anyway.
        ZMQ.close(router)
        ZMQ.close(dealer)
        # ZMQ.close(ctx)
    end
end

function broker_args()
    args = ArgParseSettings(prog="broker", autofix_names=true)  # turn "-" into "_" for arg names.

    @add_arg_table! args begin
        "--worker"
        arg_type = String
        metavar = "URL"
        default = "tcp://127.0.0.1:9467"
        help = "ZMQ DEALER address to connect to"

        "--client"
        arg_type = String
        metavar = "URL"
        default = "ipc:///tmp/chloe-client"
        help = "ZMQ ROUTER address to connect to"
    end

    parse_args(ARGS, args; as_symbols=true)

end
function broker_main()
    args = broker_args()
    Sys.set_process_title("chloe-broker")
    start_broker(args[:worker], args[:client])
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    import .Broker: broker_main
    broker_main()
end