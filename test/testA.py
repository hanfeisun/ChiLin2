execfile('../chilin_env/bin/activate_this.py', dict(__file__='../chilin_env/bin/activate_this.py'))
from chilin2.ChiLin2 import main
import sys

a = "run -c ./testA/testA.conf -v 2".split(" ") + sys.argv[1:]

print(a)
main(a)
