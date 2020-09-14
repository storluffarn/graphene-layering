
from os import system

print ("executing simulation")
#system("g++ -std=c++17 -fopenmp slip_details.cpp -o slip_details -O3 && ./slip_details")
#system("g++ -std=c++17 modtomlinjenny.cpp -o modtomlinjenny -O3 && ./modtomlinjenny")
system("g++ -std=c++17 modtomlin.cpp -o modtomlin -O3 -lgsl -lgslcblas -lm && ./modtomlin")
#system("g++ -std=c++17 modtomlinjenny.cpp -o modtomlinjenny -O3 -lgsl -lgslcblas -lm && ./modtomlinjenny")
#    
print("simulation finished")
    
print("generating plots")
    
system("python3 plotspaper2.py")
    
print("plotting finished")
#
#for i in range(100):
#    system("g++ -std=c++17 modtomlinjenny.cpp -o modtomlinjenny -O3 && ./modtomlinjenny")
#    
#    print("execution finished")
#    
#    print("generating plots")
#   
#    system("python3 plotspaper2.py")
#   
#    print("plotting finished")
#
#    print("finished {} runs".format(i))
#
