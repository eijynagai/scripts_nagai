# Useful commands in Python

# .to_string() 
# this command attached to a long dataframe allows the print command to
# shows the whole list in the screen. Usually the long list are cut and
# only the beginning and the end are shown.
print(df.to_string())


#########################
# Path() from pathlib
#########################
# Alternative to os.path
# Create a path object (PosixPath) instead of string that represent path
# This will directly use the returned object
from pathlib import Path
Path('path/file.txt').exists()
Path('new_dir').mkdir(exist_ok=True)
Path('sample_data').iterdir()

# Reading and writing with Path
f = Path('path/to/test.txt')
f.write_text('something will be included.')
f.read_text()

# Handling the metatada of a file
# for relative path:
print(f'This is the path: {f}')
# for absolute path:
print(f'This is the absolute path: {f.resolve()}')
# for the name only
print(f'File name without extension: {f.stem}')
# for the extension only
print(f'This is the extension of the file: {f.suffix}')
# for metadata stats
print(f'File stats: {f.stat()}')
# for specific metadata, such as the file size
print(f'File size (Bytes): {f.stat().st_size}')

# glob
# We can use wildcard with Path() to list all files in a directory
# The same cannot be achieved using OS library
from glob import glob
list(glob(os.path.join('directory', '*.csv'))) 


