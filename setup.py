from distutils.core import setup

setup(
    name='chilin2',
    version='0.1',
    packages=['chilin2', 'chilin2.function_template'],
    url='',
    license='MIT',
    author='Hanfei Sun, Shenglin Mei, Qian Qin ',
    author_email='',
    description='',
    scripts = ["chilin2/ChiLin2.py"],
    requires=["samflow", "jinja2"],
    package_data = {"chilin2" : ["db/ChiLinQC.db", "jinja_template/*.jinja2","db/*.txt"]}
)
