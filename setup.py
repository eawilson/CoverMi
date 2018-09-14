from setuptools import setup
import os

with open(os.path.join(os.path.dirname(__file__), "covermi", "include.py"), "rU") as f:
    exec(f.read())

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
    entry_points={"console_scripts" : ["covermi=covermi.covermigui:main",
                                       "annotatetsca=covermi.annotatetsca:main",
                                       "bringmiup=covermi.bringmiup:main",
                                      ]
                 }
)

