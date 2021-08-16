import os
import logging

try:
    # CASA 6
    logging.debug("Importing CASAtools")
    import casatools
    tb = casatools.table()
    casa6 = True

except ImportError:
    # CASA 5
    logging.debug("Import casa6 errors. Trying CASA5...")
    from taskinit import tbtool
    tb = tbtool()
    casa5 = True

class Weblog():
    def __init__(self, taskname, localdict):
        self.localdict = localdict
        self.taskname = taskname
        self.test_counter = 1
        self.all_passed = True
        #self.html = open("test_{}_weblog.html".format(self.taskname.lower()), 'w')

    #def __enter__(self):
    #    return self

    #def __exit__(self, *args):
    #    self.html.close()

    def generate_header(self, testname):
        with open("test_{}_weblog.html".format(self.taskname.lower()), 'w') as self.html:
            self.html.write('<!doctype html>' + '\n')
            self.html.write('<html lang="en">' + '\n')
            self.html.write('<head>' + '\n')
            self.html.write('<meta charset="utf-8">' + '\n')
            self.html.write('<title>{}</title>'.format(testname) + '\n')
            self.html.write('<meta name="description" content="The HTML5 Herald">' + '\n')
            self.html.write('<meta name="author" content="SitePoint">' + '\n')
            self.html.write('<link rel="stylesheet" href="css/styles.css?v=1.0">' + '\n')
            self.html.write('</head>' + '\n')
            self.html.write('<body>' + '\n')
            self.html.write('<script src="js/scripts.js"></script>' + '\n')
            self.html.write('<h1>{}</h1>'.format(testname) + '\n')

    def write_modal_style(self):
        with open("test_{}_weblog.html".format(self.taskname.lower()), 'a+') as self.html:
            self.html.write(' /* Style the Image Used to Trigger the Modal */' + '\n')
            self.html.write('.myImg {' + '\n')
            self.html.write('border-radius: 5px; cursor: pointer; transition: 0.3s; }'+ '\n')
            self.html.write('.myImg:hover{' + '\n')
            self.html.write('opacity: 0.7;}'+ '\n')
            self.html.write('/* The Modal (background) */' + '\n')
            self.html.write('.modal {' + '\n')
            self.html.write('  display: none; /* Hidden by default */' + '\n')
            self.html.write('  position: fixed; /* Stay in place */' + '\n')
            self.html.write('  z-index: 1; /* Sit on top */' + '\n')
            self.html.write('  padding-top: 100px; /* Location of the box */' + '\n')
            self.html.write('  left: 0;' + '\n')
            self.html.write('  top: 0;' + '\n')
            self.html.write('  width: 100%; /* Full width */' + '\n')
            self.html.write('  height: 100%; /* Full height */' + '\n')
            self.html.write('  overflow: auto; /* Enable scroll if needed */' + '\n')
            self.html.write('  background-color: rgb(0,0,0); /* Fallback color */' + '\n')
            self.html.write('  background-color: rgba(0,0,0,0.9); /* Black w/ opacity */' + '\n')
            self.html.write('  word-wrap: break-word; /* Wrap Text */' + '\n')
            self.html.write('}' + '\n')
            self.html.write('/* Modal Content (Image) */' + '\n')
            self.html.write('.modal-content {' + '\n')
            self.html.write('  margin: auto;' + '\n')
            self.html.write('  display: block;' + '\n')
            self.html.write('  width: 80%;' + '\n')
            self.html.write('  max-width: 700px;' + '\n')
            self.html.write('}' + '\n')
            self.html.write('/* Caption of Modal Image (Image Text) - Same Width as the Image */' + '\n')
            self.html.write('#caption {' + '\n')
            self.html.write('  margin: auto;' + '\n')
            self.html.write('  display: block;' + '\n')
            self.html.write('  width: 80%;' + '\n')
            self.html.write('  max-width: 700px;' + '\n')
            self.html.write('  text-align: center;' + '\n')
            self.html.write('  color: #ccc;' + '\n')
            self.html.write('  padding: 10px 0;' + '\n')
            self.html.write('  height: 150px;' + '\n')
            self.html.write('}' + '\n')
            self.html.write('/* Add Animation - Zoom in the Modal */' + '\n')
            self.html.write('.modal-content, #caption {' + '\n')
            self.html.write('  animation-name: zoom;' + '\n')
            self.html.write('  animation-duration: 0.6s;' + '\n')
            self.html.write('}' + '\n')
            self.html.write('@keyframes zoom {' + '\n')
            self.html.write('  from {transform:scale(0)}' + '\n')
            self.html.write('  to {transform:scale(1)}' + '\n')
            self.html.write('}' + '\n')
            self.html.write('/* The Close Button */' + '\n')
            self.html.write('.close {' + '\n')
            self.html.write('  position: absolute;' + '\n')
            self.html.write('  top: 15px;' + '\n')
            self.html.write('  right: 35px;' + '\n')
            self.html.write('  color: #f1f1f1;' + '\n')
            self.html.write('  font-size: 40px;' + '\n')
            self.html.write('  font-weight: bold;' + '\n')
            self.html.write('  transition: 0.3s;' + '\n')
            self.html.write('}' + '\n')
            self.html.write('.close:hover,' + '\n')
            self.html.write('.close:focus {' + '\n')
            self.html.write('  color: #bbb;' + '\n')
            self.html.write('  text-decoration: none;' + '\n')
            self.html.write('  cursor: pointer;' + '\n')
            self.html.write('}' + '\n')
            self.html.write('/* 100% Image Width on Smaller Screens */' + '\n')
            self.html.write('@media only screen and (max-width: 700px){' + '\n')
            self.html.write('  .modal-content {' + '\n')
            self.html.write('    width: 100%;' + '\n')
            self.html.write('  }' + '\n')
            self.html.write('} ' + '\n')



    def generate_status_table_style(self,dictionary):
        with open("test_{}_weblog.html".format(self.taskname.lower()), 'a+') as self.html:
            self.html.write('<style type="text/css">' + '\n')
            self.html.write('.collapsible {background-color: #777;color: white;cursor: pointer;padding: 18px;width: 100%;border: none;text-align: left;outline: none;font-size: 15px;}' + '\n')
            self.html.write('.active, .collapsible:hover { background-color: #555;}' + '\n')
            self.html.write('.content {padding: 0 18px;display: none;overflow: hidden;background-color: #f1f1f1;}' + '\n')
            self.html.write(".boxed { border: 1px solid green ;padding: 0 5px 0 5px;margin: 50px}" + '\n')
            self.html.write('.tg  {border-collapse:collapse;border-spacing:0;}' + '\n')
            self.html.write('.tg td{font-family:Arial, sans-serif;font-size:14px;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:black;word-wrap:break-word;}' + '\n')
            self.html.write('.tg th{font-family:Arial, sans-serif;font-size:14px;font-weight:normal;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:black;word-wrap:break-word;}' + '\n')
            self.html.write('.tg .tg-0lax{text-align:left;vertical-align:top}' + '\n')
            self.html.write('.tg .tg-ck9b{background-color:#009901;color:#32cb00;text-align:left;vertical-align:top}' + '\n')
            self.html.write('.tg .tg-r50r{background-color:#cb0000;text-align:left;vertical-align:top}' + '\n')
            self.html.write('.tg-sort-header::-moz-selection{background:0 0}.tg-sort-header::selection{background:0 0}.tg-sort-header{cursor:pointer}.tg-sort-header:after{content:'';float:right;margin-top:7px;border-width:0 5px 5px;border-style:solid;border-color:#404040 transparent;visibility:hidden}.tg-sort-header:hover:after{visibility:visible}.tg-sort-asc:after,.tg-sort-asc:hover:after,.tg-sort-desc:after{visibility:visible;opacity:.4}.tg-sort-desc:after{border-bottom:none;border-width:5px 5px 0}' + '\n')
        Weblog(self.taskname, self.localdict).write_modal_style()
        with open("test_{}_weblog.html".format(self.taskname.lower()), 'a+') as self.html:
            self.html.write('</style>' + '\n')

    def generate_tail(self,dictionary):
        with open("test_{}_weblog.html".format(self.taskname.lower()), 'a+') as self.html:
            self.html.write('<script>' + '\n')
            self.html.write('var coll = document.getElementsByClassName("collapsible");' + '\n')
            self.html.write('var i;' + '\n')
            self.html.write('for (i = 0; i < coll.length; i++) {' + '\n')
            self.html.write('  coll[i].addEventListener("click", function() {' + '\n')
            self.html.write('    this.classList.toggle("active");' + '\n')
            self.html.write('    var content = this.nextElementSibling;' + '\n')
            self.html.write('    if (content.style.display === "block") {' + '\n')
            self.html.write('      content.style.display = "none";' + '\n')
            self.html.write('    } else {' + '\n')
            self.html.write('      content.style.display = "block";' + '\n')
            self.html.write('    }' + '\n')
            self.html.write('  });' + '\n')
            self.html.write('}' + '\n')
            self.html.write('// Get the modal' + '\n')
            self.html.write("var modal = document.getElementById('myModal');" + '\n')
            self.html.write('// Get the image and insert it inside the modal - use its "alt" text as a caption' + '\n')
            self.html.write("var img = $('.myImg');" + '\n')
            self.html.write('var modalImg = $("#img01");' + '\n')
            self.html.write('var captionText = document.getElementById("caption");' + '\n')
            self.html.write("$('.myImg').click(function(){" + '\n')
            self.html.write('    modal.style.display = "block";' + '\n')
            self.html.write('    var newSrc = this.src;' + '\n')
            self.html.write("    modalImg.attr('src', newSrc);" + '\n')
            self.html.write('    captionText.innerHTML = this.alt;' + '\n')
            self.html.write('});' + '\n')
            self.html.write('// Get the <span> element that closes the modal' + '\n')
            self.html.write('var span = document.getElementsByClassName("close")[0];' + '\n')
            self.html.write('// When the user clicks on <span> (x), close the modal' + '\n')
            self.html.write('span.onclick = function() {' + '\n')
            self.html.write('  modal.style.display = "none";' + '\n')
            self.html.write('}' + '\n')
            self.html.write('</script>' + '\n')
            self.html.write("</body>" + '\n')
            self.html.write("</html>" + '\n')

    def generate_table_row(self, counter, test, description, runtime, status_color):
        with open("test_{}_weblog.html".format(self.taskname.lower()), 'a+') as self.html:
            self.html.write("<tr>" + "\n")
            self.html.write('<td class="tg-0lax">{}</td>'.format(counter) + '\n')
            self.html.write('<td class="tg-0lax">{}</td>'.format(test) + '\n')
            self.html.write('<td class="tg-0lax">{}</td>'.format(description) + '\n')
            #print("######################## RUNTIME: {}".format(runtime))
            if isinstance(runtime, float):
                self.html.write('<td class="tg-0lax">{}s</td>'.format(round(runtime,2)) + '\n')
            else:
                self.html.write('<td class="tg-0lax">{}s</td>'.format(round(runtime[0],2)) + '\n')
            self.html.write('<td class={}></td>'.format(status_color) + '\n')
            self.html.write("</tr>" + "\n")

    def generate_status_table(self, dictionary, show_passed):
        with open("test_{}_weblog.html".format(self.taskname.lower()), 'a+') as self.html:
            self.html.write('<table id="tg-cC48w" class="tg">' + '\n')
            self.html.write('<tr>' + '\n')
            self.html.write('<th class="tg-0lax">#</th>' + '\n')
            self.html.write('<th class="tg-0lax">Test Name</th>' + '\n')
            self.html.write('<th class="tg-0lax">Description </th>' + '\n')
            self.html.write('<th class="tg-0lax">Run Time</th>' + '\n')
            self.html.write('<th class="tg-0lax">Status</th>' + '\n')
            self.html.write('</tr>' + '\n')
        self.test_counter = 1
        self.all_passed = True
        for key, value in dictionary.items():
            if not show_passed:
                if dictionary[key]['status'] == True:
                    continue
            self.all_passed = False
            Weblog(self.taskname, self.localdict).generate_table_row(str(self.test_counter), str(key), dictionary[key]['description'],  dictionary[key]['runtime'], "tg-ck9b" if dictionary[key]['status'] == True else "tg-r50r" )
            self.test_counter += 1
        with open("test_{}_weblog.html".format(self.taskname.lower()), 'a+') as self.html:
            if self.all_passed:
                self.html.write("<tr>" + "\n")
                self.html.write('<td class="tg-0lax">All Test Passed</td>' + '\n')
                self.html.write("</tr>" + "\n")
            self.html.write('</table>' + '\n')


    def generate_summary_box(self, dictionary, show_passed):
        for key, value in dictionary.items():
            if not show_passed:
                if dictionary[key]['status'] == True:
                    continue
            Weblog(self.taskname, self.localdict).generate_summary_box_intro(dictionary, key, show_passed)
            Weblog(self.taskname, self.localdict).write_inline_list( dictionary[key]['taskcall'] )
            Weblog(self.taskname, self.localdict).add_miscellaneous_info(dictionary[key])
            with open("test_{}_weblog.html".format(self.taskname.lower()), 'a+') as self.html:
                self.html.write('<p><b>Re-Run:</b> {}</p>'.format(dictionary[key]['rerun'])+ '\n')
                self.html.write('</div>' + '\n')
                self.html.write('</div>' + '\n')

    def generate_summary_box_intro(self, dictionary, key, show_passed):
        with open("test_{}_weblog.html".format(self.taskname.lower()), 'a+') as self.html:

            self.html.write('<button class="collapsible">{}</button>'.format(key) + '\n')
            self.html.write('<div class="content">'+ '\n')
            self.html.write('<div class="boxed">'+ '\n')
            self.html.write('<h3>{}</h3>'.format(key) + '\n')
            self.html.write('<i><sub>{}</sub></i>'.format(dictionary[key]['description'])+ '\n')
            self.html.write('<p><b>Elapsed Time:</b> {} Seconds</p>'.format(dictionary[key]['runtime'])+ '\n')
            self.html.write('<p><b>Status:</b> {}</p>'.format("PASS" if dictionary[key]['status'] == True else "FAIL" )+ '\n')
            self.html.write('<p><b>Task Executions:</b></p>'+ '\n')

    def write_inline_list(self, array):
        with open("test_{}_weblog.html".format(self.taskname.lower()), 'a+') as self.html:
            self.html.write('<ul>' + '\n')
            for item in array:
                if isinstance(item,str):
                    self.html.write('<li>{}</li>'.format(item.replace('\\,','\n')) + '\n')
                else:
                    self.html.write('<li>{}</li>'.format(item) + '\n')
                if str(item).endswith(".png"):
                    self.html.write('<script src="https://ajax.googleapis.com/ajax/libs/jquery/2.1.1/jquery.min.js"></script>'+ '\n')
                    self.html.write('<img  class="myImg" src="{}" alt="{}" height="300" width="300">'.format(item, item) + '\n')
                    self.html.write('<div id="myModal" class="modal">'+ '\n')
                    self.html.write('  <span class="close" onclick="document.getElementById(\'myModal\').style.display=\'none\'">&times;</span>' + '\n')
                    self.html.write('  <img class="modal-content" id="img01">' + '\n')
                    self.html.write('  <div id="caption"></div>' + '\n')
                    self.html.write('</div>' + '\n')
                    self.html.write('</img>' + '\n')
            self.html.write('</ul>' + '\n')

    def write_inline_dict(self, dictionary):
        with open("test_{}_weblog.html".format(self.taskname.lower()), 'a+') as self.html:
            self.html.write('<ul>' + '\n')
            for key, value in sorted(dictionary.items()):
                self.html.write('<p><b>{}</b>: {}</p>'.format(key,value)+ '\n')
            self.html.write('</ul>' + '\n')

    def add_miscellaneous_info(self, subdictionary):

        default_keys = ['description','status','runtime','taskcall','rerun']
        for key, value in subdictionary.items():
            if key in default_keys:
                continue
            if isinstance(subdictionary[key],list):
                with open("test_{}_weblog.html".format(self.taskname.lower()), 'a+') as self.html:
                    self.html.write('<p><b>{}:</b></p>'.format(key)+ '\n')
                Weblog(self.taskname, self.localdict).write_inline_list(subdictionary[key])
            elif isinstance(subdictionary[key], dict):
                with open("test_{}_weblog.html".format(self.taskname.lower()), 'a+') as self.html:
                    self.html.write('<p><b>{}:</b></p>'.format(key)+ '\n')
                Weblog(self.taskname, self.localdict).write_inline_dict(subdictionary[key])
            else:
                with open("test_{}_weblog.html".format(self.taskname.lower()), 'a+') as self.html:
                    self.html.write('<span style="white-space: pre-line"><p><b>{}:</b> {}</p></span>'.format(key,subdictionary[key])+ '\n')


    def generate_weblog(self,show_passed):
        logging.debug("Generating Weblog: {}".format(self.taskname))
        Weblog(self.taskname, self.localdict).generate_header("Test {}".format(self.taskname))
        Weblog(self.taskname, self.localdict).generate_status_table_style(self.localdict)
        Weblog(self.taskname, self.localdict).generate_status_table(self.localdict, show_passed)
        Weblog(self.taskname, self.localdict).generate_summary_box(self.localdict, show_passed)
        Weblog(self.taskname, self.localdict).generate_tail(self.localdict)

