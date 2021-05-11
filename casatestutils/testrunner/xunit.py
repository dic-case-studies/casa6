class Xunit:
    xml_escape_table = {
        "&": "&amp;",
        '"': "&quot;",
        "'": "&apos;",
        ">": "&gt;",
        "<": "&lt;",
    }

    results = []
    fail_total = 0
    def append_result (self,testname, runtime, returncode, testerr):
        result =	{
          "testname": testname,
          "runtime" : runtime,
          "returncode": returncode,
          "testerr": testerr,
        }
        self.results.append(result)


    def xml_escape(self, text):
        return "".join(self.xml_escape_table.get(c,c) for c in text)

    def test_result_to_xml (self,result):
        self.fail_total = self.fail_total + len(result['testerr'])
        testxml = '<testcase classname="' + result['testname'].replace(".py", "") + '"' \
              + ' name="Failed tests" time="' + result['runtime'] + '">'
        if ( result['returncode'] != 0) :
            testxml = testxml + '<failure>' + str(result['testerr']) + '</failure>'
            if self.fail_total == 0:
                self.fail_total=self.fail_total + 1
        testxml = testxml + '</testcase>\n'
        return testxml

    def generateXml(self, testname):
        xmlResults = list(map(lambda result: self.test_result_to_xml(result), self.results))

        testHeader = '<?xml version="1.0" encoding="UTF-8"?>' + "\n" \
                 + '<testsuite name="UnitTests" tests="' \
                 + str(len(self.results)) + '" errors="'+ str(self.fail_total) + '"' \
                 + ' failures="' + str(self.fail_total) + '" skip="0">\n'
        print("Results: " + str(self.results))
        testFooter ="\n</testsuite>"

        # Write xUnit.xml
        xUnit = open("xUnit-"+testname+".xml", "w+")
        xUnit.write(testHeader + ''.join(xmlResults) + testFooter)
        xUnit.close()