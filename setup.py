from setuptools import setup
import os

setup(name="covermi",
    packages=["covermi"],
    version=6.0,
    description="Coverage checking for next generation sequencing panels",
    url="http://github.com/eawilson/covermi",
    author="Ed Wilson",
    author_email="edwardadrianwilson@yahoo.co.uk",
    license="MIT",
    include_package_data=True,
    zip_safe=True,
    #entry_points={"console_scripts" : ["covermi=covermi.covermimain:main",
                                       #"annotatetsca=covermi.annotatetsca:main",
                                       #"bringmiup=covermi.bringmiup:main",
                                      #]
                 #}
)

