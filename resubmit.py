import subprocess
import re
import os
import sys

def extract_ids(input_string):
    # Extract IDs from the input string
    match = re.search(r'--ids=([\d,]+)', input_string)
    if match:
        ids = match.group(1).split(',')
        return ids
    else:
        return []

def create_and_submit_condor_jobs(ids, cmd_file_path):
    try:
        with open(cmd_file_path, 'r') as cmd_file:
            cmd_content = cmd_file.read()

            # Split the content of condor.cmd by lines
            cmd_lines = cmd_content.split('\n')

            for id in ids:
                executable_path = None
                # Create a copy of the condor.cmd file
                cmd_copy_path = cmd_file_path.replace('condor.cmd', f'condor_{id}.cmd')
                with open(cmd_copy_path, 'w') as cmd_copy:
                    # Flag to check if we have removed a "queue" line
                    queue_removed = False
                    for line in cmd_lines:
                        # Modify the executable line in the copy
                        if line.strip().startswith('executable     ='):
                            executable_path = line.strip().split('=')[1].strip()
                            executable_path = executable_path.rstrip('/').replace('/condor.sh', f'/condor_{id}.sh')
                            cmd_copy.write(f"executable = {executable_path}\n")
                        elif line.strip().startswith('queue'):
                            # Remove the "queue" line
                            queue_removed = True
                        else:
                            cmd_copy.write(line + '\n')

                    # If a "queue" line was removed, add the new lines
                    if queue_removed:
                        cmd_copy.write(f"arguments = {id}\n")
                        cmd_copy.write("queue\n")
                # Print the ID of the job before submitting
                print(f"Submitting job for ID: {id}")
                # Submit the modified condor.cmd file using condor_submit -spool
                subprocess.run(["condor_submit", "-spool", cmd_copy_path])

    except FileNotFoundError:
        print(f"Error: File {cmd_file_path} not found.")

if __name__ == "__main__":
    #input_string = "bambooHTCondorResubmit --ids=254,455,456,520,521 test_12/batch/input/condor.cmd"
    if len(sys.argv) != 2:
        print("Usage: python script_name.py 'input_string'")
        sys.exit(1)

    input_string = sys.argv[1]
    ids = extract_ids(input_string)
    print("Number of jobs to be submitted: ", len(ids))
    cmd_file_path = input_string.split()[-1]

    create_and_submit_condor_jobs(ids, cmd_file_path)
