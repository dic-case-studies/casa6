import subprocess # To deploy virtual frame buffer
import tempfile
import uuid
import os
import traceback

class XvfbHelper:


    def start_virtual_frame_buffer(self):
        self.__logfile_descriptor = open("casatest_xvfb.log", 'a')
        """ Starts Xvfb and sets DISPLAY environment variable """

        def gen_displayport():
            """ produces display port in ':345' string format """
            displayport = os.getpid()
            while os.path.exists('/tmp/.X%d-lock' % displayport):
                displayport += 1

            return ":%d" % displayport

        def gen_xauth_cookie():
            """ produces cookie ready to be used with xauth add """
            try:
                cookie = subprocess.check_output(['mcookie'],
                                                 universal_newlines=True).strip()
            except Exception:
                cookie = str(uuid.uuid4()).replace('-', '')

            return cookie

        def run_xauth_xvfb(xauth_dir, xauthfile, port, cookie):
            """ runs xauth + Xvfb similarly as xvfb-run """
            subprocess.call(['xauth', '-f', xauthfile.name, 'add', port, '.', cookie],
                            stdout=self.__logfile_descriptor,
                            stderr=self.__logfile_descriptor)

            xvfb_cmd = ['Xvfb', port,'-screen','0', '2048x2048x24+32', '-auth', xauthfile.name]
            print(" ".join(xvfb_cmd))
            self.__virtual_frame_buffer_process = subprocess.Popen(
                xvfb_cmd,
                stdout=self.__logfile_descriptor, stderr=self.__logfile_descriptor,
                shell=False)

            try:
                subprocess.call(['xauth', '-f', xauthfile.name, 'remove', port],
                                stdout=self.__logfile_descriptor,
                                stderr=self.__logfile_descriptor)
                xauthfile.close()
                if os.path.isdir(xauth_dir):
                    os.rmdir(xauth_dir)
            except OSError as exc:
                print("xauth file and its subdirectory could not be removed cleanly:"
                      " {}, with exception: {}".format(xauthfile.name, exc))

        print("Starting Xvfb")
        self.__virtual_frame_buffer_port = gen_displayport()
        xauth_dir = tempfile.mkdtemp(prefix='CASA_testrunner_xauth')
        xauthfile = tempfile.NamedTemporaryFile(dir=xauth_dir)
        cookie = gen_xauth_cookie()

        try:
            run_xauth_xvfb(xauth_dir, xauthfile, self.__virtual_frame_buffer_port,
                           cookie)
            os.environ['DISPLAY'] = self.__virtual_frame_buffer_port
            print("Deployed virtual frame buffer at port {} with PID {}".format(
                         self.__virtual_frame_buffer_port,  str(self.__virtual_frame_buffer_process.pid)))
        except Exception:
            self.__virtual_frame_buffer_process = None
            formatted_traceback = traceback.format_exc()
            print("Exception deploying virtual frame buffer at {}: {}".format(
                    self.__virtual_frame_buffer_port, str(formatted_traceback)))

    def signal_stop_virtual_frame_buffer(self):
        print ("Stopping framebuffer {}".format(self.__virtual_frame_buffer_process))
        self.__logfile_descriptor.close()
        if self.__virtual_frame_buffer_process is not None:
            try:
                self.__virtual_frame_buffer_process.terminate()
                print("Virtual frame buffer deployed at {} with pid {} successfully shutdown".format(
                       self.__virtual_frame_buffer_port, self.__virtual_frame_buffer_process.pid))
            except:
                formatted_traceback = traceback.format_exc()
                print("Exception shutting down virtual frame buffer deployed at {} with pid {}: {}".format(
                        self.__virtual_frame_buffer_port,
                        str(self.__virtual_frame_buffer_process.pid),
                        str(formatted_traceback)))
        else:
            print("Virtual frame buffer (port {}) not deployed".format(self.__virtual_frame_buffer_port))

    def signal_stop_service_request(self):
        self.__stop_service_requested = True
