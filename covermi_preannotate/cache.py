from __future__ import print_function

import copy, os, pdb
import cPickle as pickle
from functools import wraps

PATH = ""
paths = {}
keyvaldata = {}

def _cachepath(kwargs):
    try:
        cachepath = kwargs["cachepath"]
        del kwargs["cachepath"]
    except KeyError:
        cachepath = None    
    return cachepath


def keyvalcache(func):
    @wraps(func)
    def wrapped(items, **kwargs):
        cachepath = _cachepath(kwargs)
        if cachepath is not None:
            funcname = func.__name__
            picklepath = os.path.join(cachepath, funcname) + ".pickle"
            if funcname not in keyvaldata:
                if os.path.exists(picklepath):
                    with open(picklepath, "rb") as f:
                        keyvaldata[funcname] = pickle.load(f)
                else:
                    keyvaldata[funcname] = {}
                    dirpath = os.path.dirname(picklepath)
                    if dirpath and not os.path.exists(dirpath):
                        os.makedirs(dirpath)
            mycache = keyvaldata[funcname]
            items = list(items)
            uncacheditems = [item for item in items if item not in mycache]

            if uncacheditems:
                for item, response in func(uncacheditems, **kwargs):
                    mycache[item] = response
                with open(picklepath, "wb") as f:
                    pickle.dump(mycache, f, -1)

            for item in items:
                try:
                    yield (item, copy.deepcopy(mycache[item]))
                except KeyError:
                    mycache[item] = None
                    yield (item, None)
        else:
            for item, response in func(items, **kwargs):
                yield (item, response)
    return wrapped


def filecache(func):
    @wraps(func)
    def wrapped(*args, **kwargs):
        cachepath = _cachepath(kwargs)
        if cachepath is None:
            raise RuntimeError
        filepath = os.path.join(cachepath, func.__name__, *args)
        dirpath = os.path.dirname(filepath)
        if not os.path.exists(filepath):
            if not os.path.exists(dirpath):
                os.makedirs(dirpath)
            func(*args, path=filepath, **kwargs)
        return filepath
    return wrapped


def objectcache(func):
    @wraps(func)
    def wrapped(*args, **kwargs):
        cachepath = _cachepath(kwargs)
        if cachepath is not None:
            filepath = os.path.join(cachepath, func.__name__, *args) + ".pickle"
            dirpath = os.path.dirname(filepath)
            try:
                with open(filepath, "rb") as f:
                    obj = pickle.load(f)
            except (IOError, EOFError):
                obj = func(*args, **kwargs)
                if not os.path.exists(dirpath):
                    os.makedirs(dirpath)
                with open(filepath, "wb") as f:
                    pickle.dump(obj, f, -1)
        else:
            obj = func(*args, **kwargs)
        return obj
    return wrapped

