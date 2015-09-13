import os, pdb





class Files(dict):

    def __init__(self, path, ext, strip_trailing_underscore=True):

        self._strip_trailing_underscore = strip_trailing_underscore

        if os.path.isfile(path):
            run, sample = os.path.split(path)
            run = os.path.basename(run)
            sample, ex = os.path.splitext(sample)
            if ex != ext:
                raise RuntimeError, "CoverMi - file "+path+" is not of type "+ext
            if self._strip_trailing_underscore:
                sample = sample.split("_")[0]
            self[(run, sample)] = path

        elif os.path.isdir(path): 
            if not hasattr(ext, "__iter__"):
                ext = (ext,)

            for root, dirnames, filenames in os.walk(path):
                break
            layout = "multipleruns"
            for filename in filenames:
                if os.path.splitext(filename)[1] in ext:
                    layout = "singlerun"
                    break
        
            if layout == "singlerun":
                self._walk_samples(path, ext[0])
            elif layout == "multipleruns":
                self._walk_runs(path, ext[0])

            for ex in ext[1:len(ext)]:
                files = type(self)(path, ex, self._strip_trailing_underscore)
                if self != files:
                    raise RuntimeError, "CoverMi - "+ext[0]+" and "+ex+" files do not match"

        else:
            raise RuntimeError, "CoverMi - "+path+" is not a file or directory"


    def _walk_runs(self, path, ext):
        for root, dirnames, filenames in os.walk(path):
            break
        for run in dirnames:
            self._walk_samples(os.path.join(root, run), ext)


    def _walk_samples(self, path, ext):
        run = os.path.basename(path)
        for root, dirnames, filenames in os.walk(path):
            for filename in filenames:
                filename, extension = os.path.splitext(filename)
                if extension == ext:
                    sample = filename.split("_")[0] if self._strip_trailing_underscore else filename
                    identifier = (run, sample)
                    if identifier in self:
                        raise RuntimeError,  "CoverMi - Sample "+sample+" duplicated in run "+run
                    self[identifier] = os.path.join(root, filename)









