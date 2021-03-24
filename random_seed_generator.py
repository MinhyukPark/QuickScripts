import random
import sys

from Bio import SeqIO
import click

@click.command()
@click.option("-n", required=True, type=int, help="number of random numbers to generate")
def random_seed_generator(n):
    """This script prints n random numbers between 1 and (2^31 - 1)
    """
    random.seed(1)
    ret_arr = [random.randint(0, (2**31) - 1) for i in range(n)]
    print(ret_arr)


if __name__ == "__main__":
    random_seed_generator()
