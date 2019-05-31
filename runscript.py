
from os import system

system("g++ -std=c++17 modtomlin.cpp -o modtomlin -O3 && ./modtomlin")

print("execution finished")

system("python tomplot.py")

print("plotting finished")

