import sys

from chilin2.ChiLin2 import main


a = "run -c ./testA/testA.conf -v 2".split(" ") + sys.argv[1:]

print(a)
main(a)
