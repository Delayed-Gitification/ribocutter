import dnaio
import pandas as pd
import argparse


def rev_c(seq):
    """
    simple function that reverse complements a given sequence
    """
    tab = str.maketrans("ACTGN", "TGACN")
    # first reverse the sequence
    seq = seq[::-1]
    # and then complement
    seq = seq.translate(tab)
    return seq


def find_guides(seq):
    gg_pos = [i for i, s in enumerate(seq) if seq[i:i + 2] == "GG" and i >= 21]  # not 22, because 0 based
    guides = []
    for i in gg_pos:
        guides.append(seq[i - 21:i - 1])
    return guides


def gen_guide_df(fastq_file, min_rl, max_rl, max_reads, max_guides):

    seqs = {}  # records which sequences are in fastq, with copy numbers
    guides_d = {}  # how many reads are targeted by each guide
    seq_guide_match = {}  # which guides match the given sequence
    total_reads = 0
    T7 = "TTCTAATACGACTCACTATA"
    overlap = "GTTTTAGAGCTAGA"

    counter = 0
    with dnaio.open(fastq_file) as f:
        for record in f:
            if min_rl < len(record.sequence) < max_rl:
                counter += 1
                try:
                    seqs[record.sequence] += 1
                except KeyError:
                    seqs[record.sequence] = 1

                total_reads += 1

                if counter >= max_reads > 0:
                    stopped_early = True
                    break

    # search for guides
    for key, value in seqs.items():
        this_guides = find_guides(key) + find_guides(rev_c(key))

        # work out how many reads each guide will target
        for guide in this_guides:
            try:
                guides_d[guide] += value
            except KeyError:
                guides_d[guide] = value

        try:
            seq_guide_match[key] += this_guides
        except KeyError:
            seq_guide_match[key] = this_guides

    # now sort this list
    sorted_guides = {k: v for k, v in sorted(guides_d.items(), key=lambda item: item[1], reverse=True)}

    # final list
    final_guides = {}
    counter = 0
    for key, value in sorted_guides.items():
        counter += 1
        if counter > max_guides:
            break
        else:
            final_guides[key] = value

    # now find how many sequences are actually targeted
    n = 0
    guide_fraction = {}
    for key, value in seqs.items():
        # find the guides associated with this seq
        this_guides = seq_guide_match[key]
        # check if they're in the list of final guides
        combined_set = set(this_guides).intersection(final_guides.keys())
        if len(combined_set) > 0:
            n += value
            for guide in combined_set:
                try:
                    guide_fraction[guide] += value / total_reads
                except KeyError:
                    guide_fraction[guide] = value / total_reads

    total_percent = 100 * n / total_reads
    print(str(total_percent) + "% of library targeted by guides")

    final_oligos = {}
    for guide in final_guides.keys():
        if not guide[0:1] == "G":
            final_oligos[guide] = T7 + "G" + guide + overlap
        else:
            final_oligos[guide] = T7 + guide + overlap

    final_df = {"oligo": list(final_oligos.values()), "target": list(final_oligos.keys()),
                "fraction": [guide_fraction[a] for a in final_oligos.keys()], "total_targeted": total_percent}

    df = pd.DataFrame.from_dict(final_df)
    return df


def main():
    # read arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", nargs='+', default=[], required=True, help="Input fastq(s)")
    parser.add_argument("-o", "--output", type=str, required=True, help="output filename")
    parser.add_argument("-r", "--max_reads", type=int, required=False, default=-1, help="max reads to examine")
    parser.add_argument("-g", "--max_guides", type=int, required=False, default=50)
    parser.add_argument("--min_read_length", type=int, required=False, default=0)
    parser.add_argument("--max_read_length", type=int, required=False, default=1000)

    args = parser.parse_args()

    if len(args.input) == 1:
        df = gen_guide_df(fastq_file=args.input[0], min_rl=args.min_read_length, max_rl=args.max_read_length,
                      max_reads=args.max_reads, max_guides=args.max_guides)
        print(df)
        df.to_csv(args.output, index=False)
    else:
        for i, filename in enumerate(args.input):
            print("Analysing " + filename)
            df = gen_guide_df(fastq_file=filename, min_rl=args.min_read_length, max_rl=args.max_read_length,
                              max_reads=args.max_reads, max_guides=args.max_guides)
            df["filename"] = filename.split("/")[-1]
            if i == 0:
                full_df = df
            else:
                full_df = full_df.append(df)
        full_df["average_fraction"] = full_df['fraction'].groupby(full_df["oligo"]).transform('sum')/len(args.input)
        print(full_df)
        full_df.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
