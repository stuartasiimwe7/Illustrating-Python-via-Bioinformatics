# %% [markdown]
# ### Preprocessing File# %%
file_name = "Insert_your_file_here.pin"

# %%
def preprocess(file_name):

    file = open(file_name, "r")
    file_write = open("new_input.pin", "w")
    lines = file.readlines()

    for a in lines:
        a = a.replace("\"", '')
        file_write.write(a)
    file_write.close()
    file.close()


# %%
preprocess(file_name) #calling the method

# %% [markdown]
# ### Calling the Percolator

# %% [markdown]
# ##### There comes a time when you need to inclide the percolator into the python code instead of running in the CMD for pre-pre-processing

# %%
import subprocess

# %%
#path to the input file
input_file = 'new_input.pin'

# %%
#command to call Percolator
call = ['percolator', input_file]

# %%
#executing the command with try & exemptions check
try:
    output = subprocess.check_output(call, universal_newlines=True)
    with open('new_output.psm','w') as file:
        file.write(output)

    print("Rescoring Successfull! File saved as \'new_output.psm\'")

except subprocess.CalledProcessError as e:
    print("RESCORING FAILED!: ",e)
