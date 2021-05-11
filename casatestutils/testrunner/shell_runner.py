import subprocess
import os

import asyncio


# Streaming logic from here:
# https://kevinmccarthy.org/2016/07/25/streaming-subprocess-stdin-and-stdout-with-asyncio-in-python/

async def _read_stream(stream, cb):
    while True:
        line = await stream.readline()
        if line:
            cb(line)
        else:
            break


async def _stream_subprocess(cmd, stdout_cb, stderr_cb, cwd):
    process = await asyncio.create_subprocess_exec(*cmd,
                                                   stdout=asyncio.subprocess.PIPE,
                                                   stderr=asyncio.subprocess.PIPE,
                                                   cwd=cwd)

    await asyncio.wait([
        _read_stream(process.stdout, stdout_cb),
        _read_stream(process.stderr, stderr_cb)
    ])
    return await process.wait()


def execute(cmd, to,stdout_cb, stderr_cb, cwd):
    loop = asyncio.get_event_loop()
    rc = loop.run_until_complete(asyncio.wait_for(
        _stream_subprocess(
            cmd,
            stdout_cb,
            stderr_cb,
            cwd
        ), timeout=to))
    #loop.close()
    return rc


class ShellRunner:
    mypath = ""


    def __init__(self):
        self.mypath = os.path.dirname(os.path.realpath(__file__))
        self.failed_tests = []
        self.test_suite_passed = False

    def process_stdout(self, x):
        message = x.decode('utf-8').rstrip("\n\r")
        print("%s" % message)

    def process_stderr(self, x):
        message = x.decode('utf-8').rstrip("\n\r")
        print("%s" % message)
        if(message.startswith("FAIL:")):
            self.failed_tests.append(message.split("FAIL: ")[1].split()[0])
        # Some MPI tests will return 1 "by design" so check for "OK" message
        if (message.startswith("OK (skipped=")):
            print("Marking test as passed.")
            self.test_suite_passed = True

    def runshell(self, cmd, timeout,cwd=None,):
        if not isinstance(cmd, list):
            raise TypeError("Input type must a list.")
        if cwd == None:
            cwd = self.mypath + "/../work"
        print("Executing: " + str(cmd))
        # result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=self.mypath + "../work")
        try: 
            result = execute(
                cmd,
                int(timeout),
                lambda x: self.process_stdout(x),
                lambda x: self.process_stderr(x),
                cwd
            )
            print ("Executor result: " + str(result))
            if (result != 0 ):
                self.failed_tests.append("Executor returned a non-zero exit code")
                if self.test_suite_passed:
                    print("Test executor returned a non-zero exit code but test script indicates a pass with" +
                        "\"OK (skipped=...\". Assuming failed MPI teardown. Marking test as passed.")
                    self.failed_tests = []
        except Exception as e:
            print(str(e))
            self.failed_tests.append("Caught exception during test execution.")
         
        print ("Failed tests:" + str(self.failed_tests))
        return self.failed_tests

