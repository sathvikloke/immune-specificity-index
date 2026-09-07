#!/usr/bin/env python3
"""Launch a command in its own session so it survives the parent shell dying.

    python3 detach.py <logfile> <cmd> [args...]

os.setsid() puts the child in a NEW session and process group, so a
process-group kill aimed at the Claude Code session cannot reach it. nohup
alone would only cover SIGHUP. Environment is inherited UNCHANGED -- in
particular no BLAS thread-count variables are set, because changing them
could change floating-point reduction order, which is the very thing the
A9 investigation is measuring.
"""
import os
import sys

log = sys.argv[1]
cmd = sys.argv[2:]

pid = os.fork()
if pid > 0:
    # Parent: report the child's pid and exit immediately.
    #
    # flush=True is load-bearing, not decoration. os._exit() skips interpreter
    # cleanup, which includes flushing stdio buffers. When stdout is a terminal
    # it is line-buffered and the pid appears anyway; when it is a PIPE -- which
    # is exactly how an agent or a script captures it -- stdout is BLOCK
    # buffered, so the pid sat unflushed in the buffer and was destroyed by
    # os._exit. Measured 2026-09-04: a detached launch printed nothing at all
    # and the pid had to be recovered from `ps`. Every caller that captured the
    # pid from this script through a pipe was silently getting an empty string.
    print(pid, flush=True)
    os._exit(0)

os.setsid()
fd = os.open(log, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
os.dup2(fd, 1)
os.dup2(fd, 2)
devnull = os.open(os.devnull, os.O_RDONLY)
os.dup2(devnull, 0)
os.environ["PYTHONUNBUFFERED"] = "1"
os.execvp(cmd[0], cmd)
