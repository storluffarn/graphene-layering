
from os import system

print("executing simulation")

#system("g++ -std=c++17 modtomlin.cpp -o modtomlin -O3 && ./modtomlin")
system("g++ -std=c++17 modtomlinjenny.cpp -o modtomlinjenny -O3 && ./modtomlinjenny")

print("execution finished")

print("generating plots")

system("python tomplot.py")

print("plotting finished")

