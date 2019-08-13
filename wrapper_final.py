import glob
import json
import os
import matplotlib.pyplot as plt
import healpy as hp
from matplotlib import cm
import matplotlib
matplotlib.use('Agg')
import ctypes
import subprocess
class CMB():
	def __init__(self):
		self.cmd = "./test.out "
		self.result_path = ""
		self.plots_path = ""
		with open('config.json') as config_file:
			data = json.load(config_file)
			self.result_path = self.result_path + str(data['output']['ouput_path'])
			self.plots_path = self.plots_path + str(data['output']['Plot_output_path'])
			#list_name = list(map(str.encode('utf-8'), list(data['input']['Files']) ))
			#list_name = data['input']['Files']
			#list_names = []
			#for i in list_name:
			#	list_names.append(str(i).encode('utf-8'))
			#print(type(data['input']['Files']))
			#print(data['input']['Files'])
			#self.arr = (ctypes.c_char_p * len(list_names))(*list_names)
			#self.cmd = self.cmd + arr
			for i in data['input']['Files']:
				self.cmd = self.cmd + ' "' + data["input"]["source_location"] + i + '"'
		print(self.cmd)
		self.mylist = [f for f in glob.glob(self.result_path + "*.fits")]

	def execute(self,input_exe = "2"):
		if str(input_exe) == "1":
			os.system(self.cmd)
		elif str(input_exe) == "2":
			if os.system('make wmap_analysis_9yr') == 0:
				os.system(self.cmd)
				#subprocess.Popen([self.cmd])
			else:
				print("Compilation error")
		else:
			print("wrong input")

	def plot(self,inp_plots="y"):
		if str(inp_plots).lower() == "y" or "yes":
			for i in self.mylist:
				print("Saving plot for :", i)
				path = self.plots_path + str(i.split('/')[-1])[:-5] + '.pdf'
				cmap = cm.bwr
				cmap = cm.binary
				cmap.set_under('blue')
				cmap.set_bad('pink')
				mapin = hp.read_map(i)
				hp.mollview(mapin, cmap=cmap)
				plt.savefig(path)
				plt.close("all")
			print("New results are saved at :", self.result_path)
			print("New plots are saved at :", self.plots_path)
		else:
			print("New results are saved at :",self.result_path)

if __name__ == "__main__":
	input_exe = input("Press 1 for old output, and 2 for new execution")
	inp_plots = input("DO you want to see maps of results (Y or N) :")
	c = CMB()
	c.execute(input_exe)
	c.plot(inp_plots)
