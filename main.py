import sys


def read_fastq(path):
    # Одно прочтение занимает четыре строки.
    with open(path) as file:
        while True:
            name = file.readline().rstrip("\r\n")
            if not name:
                break
            sequence = file.readline().rstrip("\r\n")
            plus = file.readline().rstrip("\r\n")
            quality = file.readline().rstrip("\r\n")
            yield name, sequence, plus, quality


def trim(sequence, quality, window=5, threshold=30):
    # Правило скользящего окна из Trimmomatic.
    scores = [ord(char) - 33 for char in quality]
    for i in range(len(sequence)):
        if sequence[i] == "N":
            scores[i] = 0

    if len(scores) < window:
        return "", ""
    total = sum(scores[:window])
    if total < threshold * window:
        return "", ""

    end = len(sequence)
    for i in range(len(scores) - window):
        total += scores[i + window] - scores[i]
        if total < threshold * window:
            end = i + window
            break

    # Убираем плохие основания с полученного конца.
    while end > 1 and scores[end - 1] < threshold:
        end -= 1
    return sequence[:end], quality[:end]


def write_read(file, name, sequence, plus, quality):
    file.write(f"{name}\n{sequence}\n{plus}\n{quality}\n")


def print_lengths(lengths):
    print("Минимальная длина:", min(lengths))
    print("Средняя длина:", round(sum(lengths) / len(lengths)))
    print("Максимальная длина:", max(lengths))


def main(path):
    lengths = []
    trimmed_lengths = []
    gc = 0
    quality_sum = 0
    removed = 0
    shortened = 0
    kept = 0

    with open("trimmed.fastq", "w") as trimmed_file, open("minlen60.fastq", "w") as filtered_file:
        for name, sequence, plus, quality in read_fastq(path):
            lengths.append(len(sequence))
            gc += sequence.count("G") + sequence.count("C")
            quality_sum += ord(quality[9]) - 33

            new_sequence, new_quality = trim(sequence, quality)
            if not new_sequence:
                removed += 1
                continue
            if len(new_sequence) < len(sequence):
                shortened += 1
            trimmed_lengths.append(len(new_sequence))
            write_read(trimmed_file, name, new_sequence, plus, new_quality)

            # После качества оставляем прочтения длиной от 60.
            if len(new_sequence) >= 60:
                kept += 1
                write_read(filtered_file, name, new_sequence, plus, new_quality)

    print("Всего прочтений:", len(lengths))
    print_lengths(lengths)
    print(f"GC-состав: {100 * gc / sum(lengths):.2f}%")
    print("Средний Phred в позиции 10:", round(quality_sum / len(lengths)))
    print("Полностью удалено:", removed)
    print("Укорочено и сохранено:", shortened)
    print("После тримминга:", len(trimmed_lengths))
    print_lengths(trimmed_lengths)
    print("После фильтра длины 60:", kept)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "reads.fastq")
