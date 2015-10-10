from setuptools import setup

setup(name="covermi",
    version="1.4.1",
    description="Coverage checking for next generation sequencing panels",
    url="http://github.com/eawilson/covermi",
    author="Ed Wilson",
    author_email="edwardadrianwilson@yahoo.co.uk",
    license="MIT",
    packages=["covermi"],
    include_package_data=True,
    zip_safe=True,
    entry_points={"console_scripts" : ["covermi=covermi.covermigui:main", 
                                       "covermi_make_canonical=covermi.scripts.make_canonical:main",
                                       "bringmiup=covermi.bringmiup:main",
                                       "cosmic2variants=covermi.cosmic2variants:main"]}
)

