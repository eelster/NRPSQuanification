# Hi Nathan
# Nice work! thank you for sharing!
#
# A few comments on the code below:
# a. you have zero comments in the code. at least before function definitions specify the intent and of the function. 
#    Will make it easier when someone tries to use it in the future
# b. lots of code duplication - e.g. in write_blast_standard where you basically calculate the output twice (once for s_raw and once for f).
#    better to write one function that returns a multi-line string (or list of strings) that you can later write to multiple files. 
#    that way you only build the contents of the file once.
# c. Some of the code usse elaborate string matching calculations (e.g. write_sequence). Would be simpler to use the 're' library


import os
_line_length = 80
stan_dir = "BLASTStandard"
fas_dir = "BLASTFASTA"
Fraw_dir = "BLASTFASTARaw"
Sraw_dir = "BLASTStandardRaw"
ana_dir = "Analysis"
bla_dir = "BLAST"


def simple_dir(d):
    if not os.path.exists(os.path.join(ana_dir, d)):
        os.makedirs(os.path.join(ana_dir, d))


def create_dir(f, dir):
    if not os.path.exists(os.path.join(ana_dir, os.path.join(dir, f))):
        os.makedirs(os.path.join(ana_dir, os.path.join(dir, f)))


def write_blast_standard(i, alignment, hsp, rec_dir):
    create_dir(os.path.join(stan_dir, rec_dir), bla_dir)
    create_dir(Sraw_dir, bla_dir)
    create_dir(os.path.join(Sraw_dir, rec_dir), bla_dir)
    f = open(os.path.join(ana_dir, os.path.join(bla_dir,
             os.path.join(stan_dir, os.path.join(rec_dir,
                os.path.join(alignment.title[:alignment.title.index(" ")],
                    alignment.title[:alignment.title.index(" ")] + "-" + str(i) + "-hsp"))))), 'w')
    s_raw = open(os.path.join(ana_dir, os.path.join(bla_dir,
                os.path.join(Sraw_dir, os.path.join(rec_dir,
                    alignment.title[:alignment.title.index(" ")] + "-" + str(i) + "-hsp")))), 'w')
    f.write("sequence: " + str(alignment.title[:alignment.title.index(" ")] + "-" + str(i)) + '\n')
    s_raw.write("sequence: " + str(alignment.title[:alignment.title.index(" ")] + "-" + str(i)) + '\n')
    f.write("length: " + str(abs(hsp.sbjct_end-hsp.sbjct_start)) + '\n')
    s_raw.write("length: " + str(abs(hsp.sbjct_end - hsp.sbjct_start)) + '\n')
    write_sequence('description: ', f, s_raw, alignment)
    f.write("e: " + str(hsp.expect) + '\n')
    s_raw.write("e: " + str(hsp.expect) + '\n')
    f.write("weighted e: " + str(hsp.expect * (pow(10, 37))) + '\n')
    s_raw.write("weighted e: " + str(hsp.expect * (pow(10, 37))) + '\n')
    f.write("score: " + str(hsp.score) + '\n')
    s_raw.write("score: " + str(hsp.score) + '\n')
    f.write(str(hsp.query) + '\n')
    s_raw.write(str(hsp.query) + '\n')
    f.write(str(hsp.match) + '\n')
    s_raw.write(str(hsp.match) + '\n')
    f.write(str(hsp.sbjct) + '\n')
    s_raw.write(str(hsp.sbjct) + '\n')


def write_blast_fasta(i, alignment, hsp, rec_dir):
    create_dir(os.path.join(fas_dir, rec_dir), bla_dir)
    create_dir(Fraw_dir, bla_dir)
    create_dir(os.path.join(Fraw_dir, rec_dir), bla_dir)
    f = open(os.path.join(ana_dir, os.path.join(bla_dir,
            os.path.join(fas_dir, os.path.join(rec_dir,
                os.path.join(alignment.title[:alignment.title.index(" ")],
                    alignment.title[:alignment.title.index(" ")] + "-" + str(i) + "-hsp"))))), 'w')
    f_raw = open(os.path.join(ana_dir, os.path.join(bla_dir,
                os.path.join(Fraw_dir, os.path.join(rec_dir,
                    alignment.title[:alignment.title.index(" ")] + "-" + str(i) + "-hsp")))), 'w')
    write_sequence('<', f, f_raw, alignment)
    k = 0
    while len(hsp.sbjct[k:len(hsp.sbjct)]) >= _line_length:
        f.write(hsp.sbjct[k:k+_line_length] + '\n')
        f_raw.write(hsp.sbjct[k:k + _line_length] + '\n')
        k += _line_length
    f.write(hsp.sbjct[k:len(hsp.sbjct)])
    f_raw.write(hsp.sbjct[k:len(hsp.sbjct)])


def write_full_standard(k, alignment, record):
    f = open(os.path.join(ana_dir, os.path.join(bla_dir, os.path.join(stan_dir, os.path.join(record.name, record.name + "-" + str(k) + "-alignment-" + alignment.title[:alignment.title.index(" ")])))), 'w')
    f.write("sequence: " + str(alignment.title[:alignment.title.index(" ")]) + '\n')
    f.write("length: " + str(abs(alignment.length)) + '\n')
    f.write("description: " + str(alignment.hit_def) + '\n')


def write_sequence(header, main, raw, alignment):
    if alignment.hit_def[len(alignment.hit_def) - 12:len(alignment.hit_def)] == "complete cds":
        main.write(header + str(alignment.hit_def[:len(alignment.hit_def) - 12]) + "partial cds" + '\n')
        raw.write(header + str(alignment.hit_def[:len(alignment.hit_def) - 12]) + "partial cds" + '\n')
    elif alignment.hit_def[len(alignment.hit_def) - 17:len(alignment.hit_def)] == "complete sequence":
        main.write(header + str(alignment.hit_def[:len(alignment.hit_def) - 17]) + "partial sequence" + '\n')
        raw.write(header + str(alignment.hit_def[:len(alignment.hit_def) - 17]) + "partial sequence" + '\n')
    elif alignment.hit_def[len(alignment.hit_def) - 15:len(alignment.hit_def)] == "complete genome":
        main.write(header + str(alignment.hit_def[:len(alignment.hit_def) - 15]) + "partial genome" + '\n')
        raw.write(header + str(alignment.hit_def[:len(alignment.hit_def) - 15]) + "partial genome" + '\n')
    elif alignment.hit_def[len(alignment.hit_def) - 13:len(alignment.hit_def)] == "complete mRNA":
        main.write(header + str(alignment.hit_def[:len(alignment.hit_def) - 13]) + "partial mRNA" + '\n')
        raw.write(header + str(alignment.hit_def[:len(alignment.hit_def) - 13]) + "partial mRNA" + '\n')
    else:
        main.write(header + str(alignment.hit_def) + '\n')
        raw.write(header + str(alignment.hit_def) + '\n')
