def calculate_dna(line: str) -> str:
    return f'{line.count("A")} {line.count("C")} {line.count("G")} {line.count("T")}'


if __name__ == "__main__":
    inp = input()
    print(calculate_dna(inp))
