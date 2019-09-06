
from os import system

print("executing simulation")

system("g++ -std=c++17 modtomlin.cpp -o modtomlin -O3 && ./modtomlin")

print("execution finished")

print("generating plots")

system("python tomplot.py")

print("plotting finished")

