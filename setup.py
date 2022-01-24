"""
splore

Scroll through and exploring data sets of molecules
"""
import os
import shutil
import subprocess
import sys

from setuptools import Command, find_packages, setup

import versioneer


class BuildGUICommand(Command):
    """A custom command to build the GUI and include it in the packages ``data``
    directory.
    """

    user_options = []

    def initialize_options(self) -> None:
        pass

    def finalize_options(self) -> None:
        pass

    description = "build the GUI"

    def run(self):
        cwd = os.getcwd()

        try:
            os.chdir("frontend")
            subprocess.check_call(["npm", "install"])
            subprocess.check_call(
                [
                    "npm",
                    "run",
                    "build",
                    "--",
                    "--output-path",
                    os.path.join(os.pardir, "splore", "_static"),
                    "--deploy-url",
                    "static/",
                ]
            )
            shutil.copyfile(
                os.path.join(os.pardir, "splore", "_static", "3rdpartylicenses.txt"),
                os.path.join(os.pardir, "LICENSE-3RD-PARTY"),
            )

        finally:
            os.chdir(cwd)


short_description = __doc__.split("\n")

# from https://github.com/pytest-dev/pytest-runner#conditional-requirement
needs_pytest = {"pytest", "test", "ptr"}.intersection(sys.argv)
pytest_runner = ["pytest-runner"] if needs_pytest else []

try:
    with open("README.md", "r") as handle:
        long_description = handle.read()
except IOError:
    long_description = "\n".join(short_description[2:])


setup(
    name="splore",
    author="Simon Boothroyd",
    description=short_description[0],
    long_description=long_description,
    long_description_content_type="text/markdown",
    version=versioneer.get_version(),
    cmdclass={
        **versioneer.get_cmdclass(),
        "build_gui": BuildGUICommand,
    },
    license="MIT",
    packages=find_packages(),
    include_package_data=True,
    setup_requires=[] + pytest_runner,
    entry_points={
        "console_scripts": [
            "splore=splore.cli:main",
        ],
    },
)
