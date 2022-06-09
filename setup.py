from setuptools import setup
import os

with open(os.path.join(os.path.dirname(__file__), "covermi", "version.py")) as f_in:
    exec(f_in.read())

setup(name="covermi",
    packages=["covermi"],
    version=__version__,
    description="Coverage checking for next generation sequencing panels",
    url="http://github.com/eawilson/covermi",
    author="Ed Wilson",
    author_email="edwardadrianwilson@yahoo.co.uk",
    license="MIT",
    include_package_data=True,
    zip_safe=True,
    entry_points={"console_scripts" : ["genename_bed=covermi.genename_bed:main",
                                       "design_report=covermi.design_report:main",]}
                                       #"bringmiup=covermi.bringmiup:main",
                                      #]
                 #}
)

