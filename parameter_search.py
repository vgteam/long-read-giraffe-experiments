import numpy as np
from os.path import exists
from bidict import bidict
import argparse


HASH_TO_PARAMETERS_FILE = "./hash_to_parameters.tsv"
CONFIG_FILE = "parameter_search_config.tsv"

'''
The Parameter class holds information for one giraffe parameter
It knows the range of potential values and can sample a value using different sampling functions

name is the name of the parameter getting called in giraffe (--name value)
value_range is a tuple of [range_start, range_end)
'''
class Parameter:
    def __init__(self, name, datatype, min_val, max_val, default, sampling_strategy):
        self.name = name
        self.datatype = datatype
        self.min_val = min_val
        self.max_val = max_val
        self.default = default
        self.sampling_strategy = sampling_strategy
    
    def sample(self):
        if self.min_val == self.max_val:
            return self.min_val
        elif self.datatype=="int":
            if self.sampling_strategy == "uniform":
                return np.random.randint(self.min_val, self.max_val)
            elif self.sampling_strategy == "log":
                #TODO: This is kinda hacky but it overflows if I try to exponentiate

                log_range = list(range(int(np.floor(np.log10(self.min_val))), int(np.ceil(np.log10(self.max_val)))))
                weights = [1 for x in log_range]
                weights[0] = 1- (self.min_val / np.power(10, np.ceil(np.log10(self.min_val))))
                weights[-1] = self.max_val / np.power(10, np.ceil(np.log10(self.max_val)))
                sum_weights = sum(weights)
                probs = [x / sum_weights for x in weights]

                exp_val = np.random.choice(log_range, p=probs)
                return np.random.randint(max(self.min_val, np.power(10, exp_val)), min(self.max_val, np.power(10, exp_val+1)))
                
            else: 
                print("No sampling strategy " + self.sampling_strategy + " for type " + self.datatype)
        elif self.datatype=="float":
            decimal_places = max(2, int(1/(self.max_val - self.min_val)))
            if self.sampling_strategy == "uniform":
                return round(np.random.uniform(self.min_val, self.max_val), decimal_places)
            else: 
                print("No sampling strategy " + self.sampling_strategy + " for type " + self.datatype)
    
    def __repr__(self):
        return self.name + ":\n\ttype:" + self.datatype + "\n\trange: " + str(self.min_val) + "-" + str(self.max_val) + "\n\tdefault value: " + str(self.default) + "\n\tsampling strategy: " + self.sampling_strategy
    def __str__(self):
        return self.name + ":\n\ttype:" + self.datatype + "\n\trange: " + str(self.min_val) + "-" + str(self.max_val) + "\n\tdefault value: " + str(self.default) + "\n\tsampling strategy: " + self.sampling_strategy


'''
ParameterSearch is used to store information about a set of parameters.
Define the parameters to be searched in the parameter config file (default is CONFIG_FILE)
This must be a tsv with values:
#name   type    min_val max_val default sampling_strategy
Where name is the name of the flag that giraffe uses
type is the data type (int or float)
min and max val are the range of values that the parameter can take
default is the default value from giraffe. This is used to unify old runs missing parameters
sampling_strategy is how we sample the values from the range ("uniform", "log")

Randomly sample the parameter space with sample_parameter_space(), giving it the number of sets to return.
  This will write the sampled parameters to hash_to_parameters_file
  This gets exposed to the command line with add_random_parameters.py
Load previously generated parameter sets from hash_to_parameters_file with load_parameters_from_file().
  This can also be used to run a specific parameter set. Manually write to the hash_to_parameters_file and 
  use "." as a placeholder for the hash value. It will get filled in automatically.
  This automatically gets run every time ParameterSearch is initialized.
get_hashes_and_parameter_strings() is a generator for returning a tuple of hash value and parameter string for running
  in giraffe. It returns everything stored in the ParameterSearch from sample_parameter_space() and 
  load_parameters_from_file()
get_hashes() is a generator just for the hashes
'''
class ParameterSearch:
    def __init__(self, config=CONFIG_FILE, hash_to_parameters_file = HASH_TO_PARAMETERS_FILE):

        self.hash_to_parameters_file = hash_to_parameters_file

        #This defines all parameters that can be searched and their potential values
        #It gets loaded from the config file
        self.parameters = []

        f = open(config)
        for line in f:
            if line[0] != "#":
                l = line.split()
                self.parameters.append(Parameter(l[0], l[1], 
                                                 int(l[2]) if l[1] == "int" else float(l[2]), 
                                                 int(l[3]) if l[1] == "int" else float(l[3]), 
                                                 int(l[4]) if l[1] == "int" else float(l[4]), 
                                                 l[5]) )
        f.close()

        #This maps a hash string to the set of parameters it represents, as a list of parameter values,
        # one for each parameter in self.parameters
        #Note that this is a two way dictionary. This is in case we later add parameters that didn't exist
        #before, we still want to be able to use the old hash value and fill in the default parameter values,
        #instead of re-hashing and losing the old results
        self.hash_to_parameters = bidict()

        #Initialize the parameter file, if it doesn't already exist
        if not exists(self.hash_to_parameters_file):
            f = open(self.hash_to_parameters_file, "w")
            f.write("#hash\t" + '\t'.join(param.name for param in self.parameters))
            f.close()
        else:
            #If it does exist, assume that we want it and initialize with parameters from the file
            self.load_parameters_from_file()

    #We may have previously run a parameter search and gotten hashes for the parameters run
    #If the hash file exists, load the parameter hashes from the file
    #Otherwise, start an empty file with a header for the parameters we will search
    def load_parameters_from_file(self):

        f = open(self.hash_to_parameters_file)
        
        #First, make sure that the file really holds the correct parameters
        header = f.readline().split()
        assert(header[0] == "#hash")
        for i in range(len(self.parameters)):
            assert(header[i+1] == self.parameters[i].name)
        
        rewrite_param_file = False;
        if len(header) > len(self.parameters)+1:
            #If we have defined more parameters than are in the file, then we need to re-write the 
            #file to include the new parameters with the default values
            rewrite_param_file = True
        
        #Get the values, filling in default values for things we missed
        for line in f:
            l = line.split()
            new_params = tuple([(int(l[i+1]) if self.parameters[i].datatype == "int" else float(l[i+1])) if i < len(l)-1 else self.parameters[i].default for i in range(len(self.parameters))])

            #If there wasn't a hash value for the parameter set, then make one and rewrite everything
            hash_val = l[0]
            if hash_val == ".":
                hash_val = self.parameter_tuple_to_hash(new_params) 
                rewrite_param_file = True
            self.hash_to_parameters[hash_val] = new_params
        f.close()
        
        if rewrite_param_file:
        
            f = open(self.hash_to_parameters_file, "w")
            f.write("#hash\t" + '\t'.join(param.name for param in self.parameters) + "\n")
            for k,v in self.hash_to_parameters:
                f.write(k+"\t" + '\t'.join(v))
            f.close()



    #Given a tuple representing a set of parameters, return the hash as a string
    #TODO: Idk about this...
    def parameter_tuple_to_hash(self, parameter_tuple):
        if parameter_tuple in self.hash_to_parameters.inverse:
            return self.hash_to_parameters.inverse[parameter_tuple]
        else:
            return str(abs(hash(parameter_tuple)))[:20];
        
    #Given a tuple representing a set of parameters, return a string of options to be run in giraffe
    def parameter_tuple_to_parameter_string(self, parameter_tuple):
        assert(len(parameter_tuple) == len(self.parameters))
        param_string = ""
        for i in range(len(parameter_tuple)):
            param_string+="--" + self.parameters[i].name
            param_string+=" " + str(parameter_tuple[i])
            if ( i != len(parameter_tuple)-1):
                param_string+=" "
        return param_string

    def hash_to_parameter_string(self, hash_val):
        return self.parameter_tuple_to_parameter_string(self.hash_to_parameters[hash_val])


    #Sample the parameter space and write the new parameters to HASH_TO_PARAMETERS
    def sample_parameter_space(self, count):
        f = open(self.hash_to_parameters_file, "a")
        for i in range(count):
            parameter_tuple = tuple([param.sample() for param in self.parameters])
            hash_val = self.parameter_tuple_to_hash(parameter_tuple)
            self.hash_to_parameters[hash_val] = parameter_tuple
            f.write("\n" + hash_val + "\t" + '\t'.join([str(x) for x in parameter_tuple]))
        f.close()
    
    def get_hashes(self):
        hashes = []
        for hash_val, parameter_tuple in self.hash_to_parameters.items():
            hashes.append(hash_val)
        return hashes

def main():
    parser = argparse.ArgumentParser(description="Add randomly sampled parameters to the file of parameters to search")
    parser.add_argument('--config_file', default=CONFIG_FILE, help="Config file for which parameters to sample and how") 
    parser.add_argument('--output_file', default=HASH_TO_PARAMETERS_FILE, help="File holding the parameter sets to search and their identifying hash value")
    parser.add_argument('--count', type=int, default=1000, help="How many parameters sets to sample [1000]")

    args = parser.parse_args()

    param_search = ParameterSearch(args.config_file, args.output_file)
    param_search.sample_parameter_space(args.count)


if __name__ == "__main__":
    main()
