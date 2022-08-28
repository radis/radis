import atexit
import inspect
import sys
from functools import partial
from io import StringIO


class CodeProfiler:
    def __init__(
        self,
        test_func=None,
        whitelist={"radis"},  # Look for these words in the file path.
        exclusions={"<"},  # Ignore <listcomp>, etc. in the function name.
        filename="code_profiler_out.txt",
    ):

        self.dbg_out = StringIO()
        self.test_func = test_func
        sys.setprofile(self.tracefunc)
        self.test_func = test_func
        self.memorized = set()
        self.stack_level = 0
        self.white_list = whitelist
        self.exclusions = exclusions
        atexit.register(partial(self.dump, filename=filename))

    def tracefunc(self, frame, event, arg):

        if event == "call":
            self.stack_level += 1

            unique_id = frame.f_code.co_filename + str(frame.f_lineno)
            if unique_id in self.memorized:
                return

            # Part of filename MUST be in white list.
            if any(x in frame.f_code.co_filename for x in self.white_list) and not any(
                x in frame.f_code.co_name for x in self.exclusions
            ):

                if "self" in frame.f_locals:
                    module_name = frame.f_locals["self"].__class__.__module__
                    class_name = frame.f_locals["self"].__class__.__name__
                    func_name = (
                        module_name + "." + class_name + "." + frame.f_code.co_name
                    )
                else:
                    module_name = inspect.getmodulename(frame.f_code.co_filename) + "."
                    func_name = module_name + frame.f_code.co_name

                func_name = self.stack_level * " ." + " " + func_name
                txt = "{: <60} # {}, {}\n".format(
                    func_name, frame.f_code.co_filename, frame.f_lineno
                )
                self.dbg_out.write(txt)
                self.memorized.add(unique_id)

        elif event == "return":
            if self.test_func != None:
                if self.test_func(frame):
                    self.dbg_out.write("\n\n~~~ TEST FUNCTION SUCCESS!!! ~~~\n\n")
                    self.test_func = None

            self.stack_level -= 1

    def dump(self, filename="code_profiler_out.txt"):
        ptr = self.dbg_out.tell()
        self.dbg_out.seek(0)
        with open(filename, "w") as f:
            f.write(self.dbg_out.read())
        self.dbg_out.seek(ptr)
