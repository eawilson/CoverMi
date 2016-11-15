import os


def load_conf():
    conf = {}
    try:
        with file(os.path.join(os.path.dirname(os.path.abspath(__file__)), "covermi.conf"), "rU") as f:
            for line in f:
                if "=" in line:
                    key, value = line.split("=")
                    conf[key.strip()] = value.strip()
    except IOError:
        pass
    if "panel_root" not in conf or not os.path.exists(conf["panel_root"]):
        conf["panel_root"] = os.path.dirname(os.path.expanduser("~"))
    if "bam_root" not in conf or not os.path.exists(conf["bam_root"]):
        conf["bam_root"] = os.path.dirname(os.path.expanduser("~"))
    return conf
