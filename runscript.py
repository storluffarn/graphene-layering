
from os import system

print("executing simulation")

system("g++ -std=c++17 -fopenmp slip_details.cpp -o slip_details -O3 && ./slip_details")
#system("g++ -std=c++17 -fopenmp modtomlinjenny.cpp -o modtomlinjenny -O3 && ./modtomlinjenny")
#system("g++ -std=c++17 modtomlin.cpp -o modtomlin -O3 && ./modtomlin")

print("execution finished")

print("generating plots")

system("python3 plotspaper2.py")

print("plotting finished")

