#!/usr/bin/env python

import sys
import os
import shutil

def main():
        filename = sys.argv[1]
        split_size = sys.argv[2]
        outdir = sys.argv[3]
        split_dataset(filename, split_size, outdir)
        return

def split_dataset(filename, split_size, outdir):
        print('Splitting up dataset into chunks of about ' + split_size + ' MB...')
        chunk_size = int(split_size) * 1024 * 1024
        written_bytes = 0
        basename = os.path.basename(filename)
        split_dir=1
        new_dir = os.path.join(outdir, str(split_dir))
        if os.path.exists(new_dir):        
                shutil.rmtree(new_dir)        
        os.mkdir(new_dir)
        try:
            fr = open(filename, 'r')
            split_file = os.path.join(new_dir, basename)
            fw = open(split_file, 'w')
            for line in fr:
                if line.startswith('>') and written_bytes >= chunk_size:
                    fw.close()
                    written_bytes = 0
                    split_dir +=1
                    new_dir = os.path.join(outdir, str(split_dir))
                    if os.path.exists(new_dir):
                            shutil.rmtree(new_dir)
                    os.mkdir(new_dir)
                    split_file = os.path.join(new_dir, basename)                    
                    fw = open(split_file, 'w')
                fw.write(line)
                written_bytes += len(line)
            fr.close()
            fw.close()
        except Exception as e:
            error = 'Failed to split the dataset! - ERROR: ' + str(e)
            print(error)
            return False
        print('Created ' + str(split_dir) + 'splits.')
        return True

if __name__ == "__main__":
        main()
