import random

import click
import numpy as np


@click.command()
@click.option("--input-filename", required=True, type=click.Path(exists=True), help="Input iqtree model file to be parsed")
@click.option("--label", required=True, type=str, help="Identifiable label for the feature vector")
@click.option("--output-filename", required=True, type=click.Path(exists=False), help="Output file to append the rates to")
@click.option("--iqtree/--no-iqtree", required=False, default=True, help="Whether to output using RAxML-ng style or IQTree style")
def extract_feature_vector(input_filename, label, output_filename, iqtree):
    ''' This program parses an IQTree file and then prints the model string
    It prints a string in a Model{0.0, 1.0....}+R{w1,r1, w2,r2, ...., wk,rk} format
    '''
    model_rates = []
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
            cur_line = f.readline()
    # print(model_family)
    # print("model_rates: ", model_rates)
    # print(model_rate_heterogeneity)
    # print(state_freq)
    # print(pinv)
    # print(alpha)
    # print(free_rates)
    # print(final_string, end="")
    with open(output_filename, "a+") as f:
        f.write(label)
        f.write(":")
        for model_rate in model_rates:
            f.write(" ")
            f.write(model_rate)
        f.write("\n")

if __name__ == "__main__":
    extract_feature_vector();


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
