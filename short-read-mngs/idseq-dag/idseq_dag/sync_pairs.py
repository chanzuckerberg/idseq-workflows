import re 
RE_SPLIT = re.compile('(/[12])?\t')

def sync_pairs(fastq_files, max_discrepant_fraction=0):
    """The given fastq_files contain the same read IDs but in different order.
    Output the same data in synchronized order. Omit any reads missing their mate
    up to max_discrepant_fraction if necessary. If more must be suppressed,
    indicate it in the second value of the returned tuple.
    """
    if len(fastq_files) != 2:
        return fastq_files, False

    output_fnames = [ifn + ".synchronized_pairs.fq" for ifn in fastq_files]
    with open(fastq_files[0], "rb") as if_0, open(fastq_files[1],
                                                    "rb") as if_1:
        with open(output_fnames[0], "wb") as of_0, open(
                output_fnames[1], "wb") as of_1:
            outstanding_r0, outstanding_r1, max_mem, total = sync_pairs_work(
                of_0, of_1, if_0, if_1)
    if max_mem:
        # This will be printed if some pairs were out of order.
        msg = "WARNING: Pair order out of sync in {fqf}. Synchronized using RAM for {max_mem} pairs."
        msg = msg.format(fqf=fastq_files, max_mem=max_mem)
        print(msg)

    discrepancies_count = len(outstanding_r0) + len(outstanding_r1)
    if discrepancies_count:
        msg = "WARNING: Found {dc} broken pairs in {fqf}, e.g., {example}."
        msg = msg.format(dc=discrepancies_count,
                            fqf=fastq_files,
                            example=(outstanding_r0 or outstanding_r1).popitem()[0])
        print(msg)
    too_discrepant = (discrepancies_count > max_discrepant_fraction * total)
    return output_fnames, too_discrepant

def extract_rid(s):
    return RE_SPLIT.split(s, 1)[0].strip()

def sync_pairs_work(of0, of1, if0, if1):
    # TODO: Use this as a template for merging fasta?
    outstanding_r0 = {}
    outstanding_r1 = {}
    mem = 0
    max_mem = 0
    total = 0
    while True:
        r0, r0id = get_read(if0)
        r1, r1id = get_read(if1)
        if not r0 and not r1:
            break
        total += 2
        if r0id == r1id:
            # If the input pairs are already synchronized, we take this
            # branch on every iteration.
            write_lines(of0, r0)
            write_lines(of1, r1)
        else:
            mem, max_mem = handle_outstanding_read(
                r0, r0id, outstanding_r0, outstanding_r1, of0, of1, mem,
                max_mem)
            mem, max_mem = handle_outstanding_read(
                r1, r1id, outstanding_r1, outstanding_r0, of1, of0, mem,
                max_mem)
    return outstanding_r0, outstanding_r1, max_mem, total

def handle_outstanding_read(r0, r0id, outstanding_r0, outstanding_r1, of0,
                            of1, mem, max_mem):
    # If read r0 completes an outstanding r1, output the pair (r0, r1).
    # Else r0 becomes outstanding, so in future some r1 may complete it.
    if r0id:
        if r0id in outstanding_r1:
            write_lines(of0, r0)
            write_lines(of1, outstanding_r1.pop(r0id))
            mem -= 1
        else:
            outstanding_r0[r0id] = r0
            mem += 1
            if mem > max_mem:
                max_mem = mem
    return mem, max_mem

def get_read(f):
    # The FASTQ/FASTA format specifies that each read consists of 4/2 lines,
    # the first of which begins with @/> followed by read ID.
    read, rid = [], None
    line = f.readline()
    if line:
        if line[0] == 64:  # Equivalent to '@', fastq format
            rid = extract_rid(line.decode('utf-8'))
            read.append(line)
            for _ in range(3):
                read.append(f.readline())
        elif line[0] == 62:  # Equivalent to '>', fasta format
            rid = extract_rid(line.decode('utf-8'))
            read.append(line)
            read.append(f.readline())
        else:
            raise RuntimeError("sync pair failed. unknown ")
    return read, rid

def write_lines(of, lines):
    for l in lines:
        of.write(l)

def unmapped_files_in(folder, num_inputs):
    return [f"{folder}/Unmapped.out.mate{i+1}" for i in range(num_inputs)]


def main():
    import sys
    pair1 = sys.argv[1]
    pair2 = sys.argv[2]
    output_files, too_discrepant = sync_pairs((pair1, pair2))
    if too_discrepant:
        raise ValueError("pairs are too discrepant")

if __name__ == "__main__":
    main()