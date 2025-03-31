import os
import re
import subprocess
import time

trace_base_path = "/mnt/storage/traces/gtrace_v2/"

# List all folders in the trace_base_path
trace_folders = [f for f in os.listdir(trace_base_path) if os.path.isdir(os.path.join(trace_base_path, f))]
print("Trace folders:")
print(trace_folders)

trace_folder_full_path = [os.path.join(trace_base_path, f) for f in trace_folders]
print("Trace folders full path:")
print(trace_folder_full_path)

core_map = {}

# Build core_map by reading the *_info.textproto file for each trace
for trace in trace_folders:
    description_file = os.path.join(trace_base_path, trace + "_info.textproto")
    with open(description_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            pattern = r"peak_live_core_count: (\d+)"
            match = re.search(pattern, line)
            if match:
                core_map[trace] = match.group(1)
                break

print("Core map:")
print(core_map)

# Ensure the output directory exists
output_dir = "./build/batch"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Create individual directories for each trace inside the output directory
for trace in trace_folders:
    trace_output_dir = os.path.join(output_dir, trace)
    os.makedirs(trace_output_dir, exist_ok=True)
    print(f"Created directory for {trace} in {output_dir}")

# Launch processes asynchronously and keep track of them along with their log file handles
processes = []  # will store tuples of (trace, process, log_file)

for trace in trace_folders:
    trace_output_dir = os.path.join(output_dir, trace)
    # Build the command
    command = [
        "./build2/test",
        "-t", os.path.join(trace_base_path, trace),
        "-c", "100000000000",
        "-n", core_map.get(trace, "1"),  # defaulting to "1" if not found
        "-q",
        "-p", trace_output_dir,
        "-f", trace
    ]
    print("Command:")
    print(command)
    
    # Set up the log file path and open the file outside the context manager so it remains open
    log_file_path = os.path.join(trace_output_dir, f"{trace}.log")
    # if os.path.exists(log_file_path):
    #     print(f"Log file {log_file_path} already exists. skipping.")
    #     continue
    
    log_file = open(log_file_path, "w")
    
    
    
    # Launch the process asynchronously, redirecting stdout and stderr to the log file
    process = subprocess.Popen(command, stdout=log_file, stderr=subprocess.STDOUT)
    processes.append((trace, process, log_file))

# Monitor processes and print a message as each one finishes
while processes:
    # Iterate over a copy of the list to safely remove items
    for trace, process, log_file in processes.copy():
        retcode = process.poll()
        if retcode is not None:  # Process has finished
            print(f"Process for trace '{trace}' finished with exit code {retcode}.")
            log_file.close()  # Close the log file for this process
            processes.remove((trace, process, log_file))
    time.sleep(0.5)  # Wait briefly before checking again

print("All processes have finished.")
