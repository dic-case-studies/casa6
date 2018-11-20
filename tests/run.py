###
### use like:
###
###      PYTHONPATH=build/lib.macosx-10.12-x86_64-3.6 python tests/run.py tests/tools/image/test_ia_imageconcat.py tests/tools/image/test_ia_makecomplex.py
###
import os
import sys
import time
import inspect
import unittest
import subprocess
from subprocess import Popen, PIPE
from distutils.dir_util import remove_tree
from xml.sax.saxutils import escape

###
### xUnit generation
###

xml_escape_table = {
    "&": "&amp;",
    '"': "&quot;",
    "'": "&apos;",
    ">": "&gt;",
    "<": "&lt;",
}

def xml_escape(text):
    return "".join(xml_escape_table.get(c,c) for c in text)

def readFile(f):
    with open(f, 'r') as f1:
      return f1.read()

def test_result_to_xml (result):

    returncode, testname, run_time, teststdout, testerr = result
    testxml = '<testcase classname="' + testname + '"' \
          + ' name="full log" time="' + str(round(run_time)) + '">'
    if ( returncode != 0) :
       testxml = testxml + '<failure>' + xml_escape(readFile(testerr)) + '</failure>'
    testxml = testxml + '</testcase>\n'
    return testxml

###
### support routines
###
def mkpath(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

pyversion = float(sys.version_info[0]) + float(sys.version_info[1]) / 10.0

if pyversion < 3:
    str_encode = str
    str_decode = str
    def pipe_decode(output):
        return output
else:
    def str_encode(s):
        return bytes(s,sys.getdefaultencoding())
    def str_decode(bs):
        return bs.decode(sys.getdefaultencoding(),"strict")
    def pipe_decode(output):
        if isinstance(output,bytes) or isinstance(output,bytearray):
            return str_decode(output)
        elif isinstance(output,tuple):
            return (str_decode(output[0]),str_decode(output[1]))
        else:
            return ("","")

test_env = os.environ.copy( )
def dump_output(working_dir, bname, out, err):
    stdout_path = "%s/%s-stdout.txt" % (working_dir,bname)
    stdout_fd = open(stdout_path, 'w')
    stdout_fd.write(out)
    stdout_fd.close( )
    stderr_path = "%s/%s-stderr.txt" % (working_dir,bname)
    stderr_fd = open(stderr_path, 'w')
    stderr_fd.write(err)
    stderr_fd.close( )
    return (stdout_path,stderr_path)

def run_test(tabwidth,test_path,working_dir):
    label = '.'.join(os.path.basename(test_path).split('.')[:-1])
    sys.stdout.write(label + '.' * (tabwidth - len(label)))
    sys.stdout.flush( )
    start_time = time.time()
    proc = Popen( [sys.executable,test_path], cwd=working_dir, env=test_env,
                  stdout=subprocess.PIPE, stderr=subprocess.PIPE )
    (output, error) = pipe_decode(proc.communicate( ))
    exit_code = proc.wait( )
    (stdout_path,stderr_path) = dump_output(working_dir,"log",output,error)
    end_time = time.time()
    run_time = end_time-start_time
    print(" ok" if exit_code == 0 else " fail")
    return (exit_code, label, run_time, stdout_path, stderr_path)

###
### collect paths and test module names
###

if len(sys.argv) == 1:
    tests = [ ]
    run_py = os.path.realpath(__file__)
    test_dir = os.path.dirname(run_py)
    for dir, subdirs, files in os.walk(test_dir):
        for f in files:
            if f.endswith(".py") and f.startswith("test_"):
                workingdir = "%s/work/%s" % (test_dir,f[:-3])
                if os.path.exists(workingdir):
                    remove_tree(workingdir)
                mkpath(workingdir)
                tests.append((os.path.abspath("%s/%s" % (dir,f)),workingdir))

    testwidth = 0 if len(tests) == 0 else max(map(lambda x: len(os.path.basename(x[0]))+3,tests))
    tabwidth = max(testwidth,45)

    start_time = time.time()

    results = list(map(lambda params: run_test(tabwidth,*params),tests))

    print('-' * (tabwidth + 8))
    passed = list(filter(lambda v: v[0] == 0,results))
    failed = list(filter(lambda v: v[0] != 0,results))

    # Construct xUnit.xml
    testHeader = '<?xml version="1.0" encoding="UTF-8"?>' + "\n" \
             + '<testsuite name="UnitTests" tests="' \
             + str(len(results)) + '" errors="0"' \
             + ' failures="' + str(len(failed)) + '" skip="0">\n'
    xmlResults = list(map(lambda result: test_result_to_xml (result), results))
    testFooter ="\n</testsuite>"

    # Write xUnit.xml
    xUnit = open("xUnit.xml","w+")
    xUnit.write(testHeader + ''.join(xmlResults) + testFooter)
    xUnit.close()

    end_time = time.time()

    print("ran %s tests in %.02f minutes, %d passed, %d failed" % (len(results),(end_time-start_time) / 60.0,len(passed),len(failed)))
    print("OK" if len(failed) == 0 else "FAIL")


    sys.exit(0 if len(failed) == 0 else 1)

else:
    test_paths = set( )
    test_modules = [ ]

    for i in sys.argv[1:]:
        ## note directory
        suite_path = os.path.dirname(i)
        test_paths.add(suite_path)
        ## discover module name
        suite_module, xxx = os.path.splitext(os.path.basename(i))
        test_modules.append(suite_module)

    sys.path = list(test_paths) + sys.path

    suite = unittest.defaultTestLoader.loadTestsFromNames(test_modules)
    unittest.TextTestRunner().run(suite)
