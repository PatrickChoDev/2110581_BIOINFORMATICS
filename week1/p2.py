def replace_rna(line: str) -> str:
    return line.replace("T","U")

if __name__ == "__main__":
    print(replace_rna(input()))
