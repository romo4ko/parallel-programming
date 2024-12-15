import matplotlib.pyplot as plt
import re

def extract_data(file_path):
    data = {}
    with open(file_path, 'r') as file:
        for line in file:
            match = re.match(r"(Sequential|Parallel) (sum|average|min/max|euclidean norm|manhattan norm|dot product): (\d+)", line)
            if match:
                mode = match.group(1)
                operation = match.group(2)
                duration = int(match.group(3))
                if operation not in data:
                    data[operation] = {'Sequential': [], 'Parallel': []}
                data[operation][mode].append(duration)
    return data

def visualize_data(data):
    for operation, timings in data.items():
        iterations = range(1, len(timings['Sequential']) + 1)
        plt.plot(iterations, timings['Sequential'], marker='o', label='Sequential')
        plt.plot(iterations, timings['Parallel'], marker='x', label='Parallel')
        plt.xlabel('Iteration')
        plt.ylabel('Time (ms)')
        plt.title(f'{operation} Performance')
        plt.legend()
        plt.grid(True)
        plt.savefig(f'{operation}_plot.png')  # Save the plot to a file
        plt.show()

if __name__ == "__main__":
    data = extract_data("results.txt")
    visualize_data(data)