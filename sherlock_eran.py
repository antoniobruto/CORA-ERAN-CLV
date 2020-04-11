import re
import numpy as np
import argparse
import os


path = "./sherlock_wt/"
for files_ in os.listdir(path):
	print(files_)
	file_name = str(files_)

	fp = open(path+file_name,"r")
	fp_len = open(path+file_name, "r")
	
	file_length = sum(1 for line in fp_len)
	print(file_length)
	
	data=""

	cnt = 1

	inputs = int(fp.readline().strip())
	outputs = int(fp.readline().strip())
	no_hidden_layers = int(fp.readline().strip())

	# print("ip : {}, op: {}, no_hidden_layers : {}".format(inputs, outputs, no_hidden_layers))

	no_neurons_layers = []


	for i in range(no_hidden_layers):
		no_neurons_layers.append(int(fp.readline().strip()))


	layers = []
	bias = []
	prev_layer_len = inputs


	for i in range(len(no_neurons_layers)):
		weights = []
		bias = []
		# print("Layer: {}, no of neurons :{}, prev_layer_len :{}".format(i, no_neurons_layers[i], prev_layer_len))
		for j in range(no_neurons_layers[i]):
			wt = []
			for k in range(prev_layer_len):
				wt.append(float(fp.readline().strip()))
			
			weights.append(wt)
			bias.append(float(fp.readline().strip()))
		
		data += "ReLU\n"
		data += str(weights)
		data += "\n"
		data += str(bias)
		data += "\n"
		# print("Layer : {}".format(prev_layer_len))
		prev_layer_len = no_neurons_layers[i]


	weights = []
	for j in range(outputs):
		wt = []
		b = []
		for k in range(prev_layer_len):
			wt.append(float(fp.readline().strip()))
		b.append(float(fp.readline().strip()))
		weights.append(wt)
		

		data += "ReLU\n"
		data += str(weights) + "\n" 
		data += str(b) + "\n"

	fp.close()

	fp_op = open('./eran_wt/'+file_name+'.tf','w')
	fp_op.write(data)
	fp_op.close()

