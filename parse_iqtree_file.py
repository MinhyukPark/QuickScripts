import random

import click
import numpy as np

MODEL_MAP = {
    "JC": "000000",
    "JC69": "000000",
    "F81": "000000",
    "K80": "010010",
    "K2P": "010010",
    "HKY": "010010",
    "TN": "010020",
    "TN93": "010020",
    "TNe": "010020",
    "K81": "012210",
    "K3P": "012210",
    "K81u": "012210",
    "TPM2": "010212",
    "TPM2u": "010212",
    "TPM3": "012012",
    "TPM3u": "012012",
    "TIM": "012230",
    "TIMe": "012230",
    "TIM2": "010232",
    "TIM2e": "010232",
    "TIM3": "012032",
    "TIM3e": "012032",
    "TVM": "012314",
    "TVMe": "012314",
    "SYM": "012345",
    "GTR": "012345",
}

@click.command()
@click.option("--input-filename", required=True, type=click.Path(exists=True), help="Input iqtree model file to be parsed")
@click.option("--iqtree/--no-iqtree", required=False, default=True, help="Whether to output using RAxML-ng style or IQTree style")
def main_entry(input_filename, iqtree):
    ''' This program parses an IQTree file and then prints the model string
    It prints a string in a Model{0.0, 1.0....}+R{w1,r1, w2,r2, ...., wk,rk} format
    '''
    model_family = ""
    model_rates = []
    model_rate_heterogeneity = []
    state_freq = []
    free_rates = []
    pinv = None
    alpha = None
    with open(input_filename) as f:
        cur_line = f.readline()
        while cur_line:
            if("Model of substitution" in cur_line):
                cur_arr = cur_line.split()
                cur_model = cur_arr[3].strip()
                model_family = cur_model.split("+")[0]
                model_rate_heterogeneity = cur_model.split("+")[1:]
                for i in range(3):
                    cur_line = f.readline()
                for i in range(6):
                    cur_line = f.readline()
                    cur_rate = cur_line.split()[1]
                    model_rates.append(cur_rate.strip())
            if("F" in model_rate_heterogeneity):
                if("State frequencies" in cur_line):
                    cur_line = f.readline()
                    for i in range(4):
                        cur_line = f.readline()
                        state_freq.append(cur_line.split("=")[1].strip())
            if any("I" in cur for cur in model_rate_heterogeneity):
                if("Proportion of invariable sites" in cur_line):
                    pinv = cur_line.split(":")[1].strip()
            if any("G" in cur for cur in model_rate_heterogeneity):
                if("Gamma shape alpha" in cur_line):
                    alpha = cur_line.split(":")[1].strip()
            if any("R" in cur for cur in model_rate_heterogeneity):
                num_cat = 0
                for cur in model_rate_heterogeneity:
                    if("R" in cur):
                        if(cur == "R"):
                            num_cat = 4
                        else:
                            num_cat = int(cur[1:])
                if("Model of rate heterogeneity" in cur_line and "FreeRate" in  cur_line):
                    cur_line = f.readline()
                    while("Category  Relative_rate  Proportion" not in cur_line):
                        cur_line = f.readline()
                    for i in range(num_cat):
                        cur_line = f.readline()
                        cur_weight_rate = cur_line.split()
                        free_rates.append(cur_weight_rate[2].strip())
                        free_rates.append(cur_weight_rate[1].strip())
                    free_weights = []
                    for free_rate in free_rates[::2]:
                        free_weights.append(float(free_rate))
                    weight_sum = sum(free_weights)
                    if(weight_sum != 1):
                        # I shouldn't have to do this correction but iqtree returns proportions that don't add up to 1 sometimes
                        normalized_weights = np.array(free_weights) / weight_sum
                        for weight_index,normalized_weight in enumerate(normalized_weights):
                            free_rates[weight_index * 2] = str(normalized_weight)
            cur_line = f.readline()
    # print(model_family)
    # print("model_rates: ", model_rates)
    # print(model_rate_heterogeneity)
    # print(state_freq)
    # print(pinv)
    # print(alpha)
    # print(free_rates)
    current_model_dict = {}
    max_position = int(MODEL_MAP[model_family][0])
    for position_index, position in enumerate(MODEL_MAP[model_family]):
        max_position = max(max_position, int(position))
        if int(position) not in current_model_dict:
            current_model_dict[int(position)] = model_rates[position_index]
    final_string = model_family
    max_position += 1
    skip_position = int(MODEL_MAP[model_family][5])
    # print("skip: ", skip_position)
    # print("model_dict: ", current_model_dict)
    separator = ","
    if not iqtree:
        separator = "/"

    if(len(current_model_dict) > 1):
        final_string += "{"
        for i in range(max_position):
            if(i != skip_position):
                final_string += (current_model_dict[i] + separator)
        final_string = final_string[:-1] + "}"

    if(len(model_rate_heterogeneity) > 0):
        if any("F" in current_rate for current_rate in model_rate_heterogeneity):
            if iqtree:
                final_string += "+F{"
            else:
                final_string += "+FU{"
            for i in range(4):
                final_string += (state_freq[i] + separator)
            final_string = final_string[:-1] + "}"
        else:
            pass
            # final_string += "+FQ"
        if any("I" in current_rate for current_rate in model_rate_heterogeneity):
            if iqtree:
                final_string += ("+I{" + pinv + "}")
            else:
                final_string += ("+IU{" + pinv + "}")
        if any("G" in current_rate for current_rate in model_rate_heterogeneity):
            current_rate_heterogeneity = ""
            for current_rate in model_rate_heterogoneity:
                if("G" in current_rate):
                    current_rate_heterogeneity = current_rate
            final_string += ("+" + current_rateheterogeneity + "{" + alpha + "}")
        if any("R" in current_rate for current_rate in model_rate_heterogeneity):
            for current_rate in model_rate_heterogeneity:
                if("R" in current_rate):
                    current_rate_heterogeneity = current_rate
            final_string += "+" + current_rate_heterogeneity + "{"
            if iqtree:
                for i in range(len(free_rates)):
                    final_string += (free_rates[i] + separator)
                final_string = final_string[:-1] + "}"
            else:
                for i in range(1, len(free_rates), 2):
                    final_string += (free_rates[i] + separator)
                final_string = final_string[:-1] + "}"
                final_string += "{"
                for i in range(0, len(free_rates), 2):
                    final_string += (free_rates[i] + separator)
                final_string = final_string[:-1] + "}"
    print(final_string, end="")

if __name__ == "__main__":
    main_entry();


'''
SYM{1.2806,3.6117,0.5469,1.1950,4.6534}+R{0.1227,0.05953,0.2267,0.2769,0.2419,0.6721,0.2795,1.367,0.1291,2.983}

TIM2{}

'''

'''
relevant file section
Model of substitution: SYM+R5

Rate parameter R:

  A-C: 1.2806
  A-G: 3.6117
  A-T: 0.5469
  C-G: 1.1950
  C-T: 4.6534
  G-T: 1.0000

State frequencies: (equal frequencies)

Rate matrix Q:

  A   -0.8853    0.2084    0.5879   0.08902
  C    0.2084     -1.16    0.1945    0.7574
  G    0.5879    0.1945   -0.9451    0.1628
  T   0.08902    0.7574    0.1628    -1.009

Model of rate heterogeneity: FreeRate with 5 categories
Site proportion and rates:  (0.1227,0.05953) (0.2267,0.2769) (0.2419,0.6721) (0.2795,1.367) (0.1291,2.983)

 Category  Relative_rate  Proportion
  1         0.05953        0.1227
  2         0.2769         0.2267
  3         0.6721         0.2419
  4         1.367          0.2795
  5         2.983          0.1291

'''
'''
SUBSTITUTION PROCESS
--------------------

Model of substitution: TVM+F+R6

Rate parameter R:

  A-C: 1.0310
  A-G: 1.3248
  A-T: 1.6824
  C-G: 0.6068
  C-T: 1.3248
  G-T: 1.0000

State frequencies: (empirical counts from alignment)

  pi(A) = 0.2034
  pi(C) = 0.2875
  pi(G) = 0.3233
  pi(T) = 0.1859

Rate matrix Q:

  A    -1.302    0.3721    0.5377    0.3926
  C    0.2632   -0.8187    0.2463    0.3092
  G    0.3383     0.219   -0.7906    0.2334
  T    0.4296    0.4782    0.4059    -1.314

Model of rate heterogeneity: FreeRate with 6 categories
Site proportion and rates:  (0.1847,0.06256) (0.1348,0.2105) (0.1586,0.4921) (0.1908,0.9559) (0.2132,1.681) (0.1179,2.895)

 Category  Relative_rate  Proportion
  1         0.06256        0.1847
  2         0.2105         0.1348
  3         0.4921         0.1586
  4         0.9559         0.1908
  5         1.681          0.2132
  6         2.895          0.1179
'''

'''
Model of substitution: SYM+I+G4

Rate parameter R:

  A-C: 1.1842
  A-G: 3.3784
  A-T: 0.5149
  C-G: 1.1312
  C-T: 4.2375
  G-T: 1.0000

State frequencies: (equal frequencies)

Rate matrix Q:

  A   -0.8872    0.2069    0.5903   0.08996
  C    0.2069    -1.145    0.1977    0.7404
  G    0.5903    0.1977   -0.9627    0.1747
  T   0.08996    0.7404    0.1747    -1.005

Model of rate heterogeneity: Invar+Gamma with 4 categories
Proportion of invariable sites: 0.01727
Gamma shape alpha: 0.9563

 Category  Relative_rate  Proportion
  0         0              0.01727
  1         0.1304         0.2457
  2         0.4714         0.2457
  3         1.01           0.2457
  4         2.459          0.2457
Relative rates are computed as MEAN of the portion of the Gamma distribution falling in the category.

MAXIMUM LIKELIHOOD TREE
'''
