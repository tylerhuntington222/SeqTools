
def write_decoding_to_file(step_data, outfile):
    with open(outfile, 'w') as f:
        f.write("# Viterbi_decoding posterior_decoding posterior_mean")
        f.write("\n")
        for vit_dec, post_dec, post_mean in step_data:
            line = " ".join([vit_dec, post_dec, post_mean])
            f.write(line)
            f.write("\n")


if __name__ == "__main__":
    test_data = []
    for i in range(10):
        test_data.append(("a", "b", "c"))
    write_decoding_to_file(test_data, "test.txt")



def write_log_liklihoods(outfile, likelihoods):
    with open(outfile, 'w') as f:
        line_1 = "# Likelihood under {initial, estimated} parameters"
        f.write(line_1)
        f.write(str(likelihoods[0]))
        f.write(str(likelihoods[1]))

def write_estimated_params_file(outfile, params):
    string = display_params(params)
    with open(outfile, 'w') as f:
        f.write(string)
