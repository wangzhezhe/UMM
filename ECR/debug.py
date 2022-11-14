# get the data and draw the plot out

import sys
import statistics 
import matplotlib.pyplot as plt


if __name__ == "__main__":
    
    if len(sys.argv)!=2:
        print("<binary> <outputname>")
        exit()

    output_file = sys.argv[1]

    # label the boundry
    fig, ax = plt.subplots(figsize=(6, 6))

    #plt.xlim([0, 6])
    #plt.ylim([0, 6])

    for line in sys.stdin:
        input_list=line.split()
        offset=3

        # get the bound for one cell
        xmin = input_list[offset+0]
        xmax = input_list[offset+1]
        ymin = input_list[offset+2]
        ymax = input_list[offset+3]

        print(xmin, ymin)
        print(xmax, ymin)
        print(xmax, ymax)
        print(xmin, ymax)

        plt.scatter(xmin, ymin)
        plt.scatter(xmax, ymin)
        plt.scatter(xmax, ymax)
        plt.scatter(xmin, ymax)

    #new_list = [(elem) for elem in line.split()]
    #if(len(new_list)>0):
    #  data_list.append(new_list[0])

#plt.scatter(0, 0)
#plt.scatter(1.5,0)

#plt.scatter(1.5, 1.5)
#plt.scatter(0, 1.5)

    plt.savefig(output_file+'png')