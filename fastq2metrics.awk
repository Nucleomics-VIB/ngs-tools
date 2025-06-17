NOT OK - TO BE CHECKED

#!/bin/awk -f

# Function to calculate average quality
#function aveQual(quals, len) {
#    sum = 0
#    for (i = 1; i <= len; i++) {
#        q = ord(substr(quals, i, 1))
#        sum += 10^(q / -10)
#    }
#    avg_prob = sum / len
#    avg_phred = -10 * log(avg_prob) / log(10)
#    return avg_phred
#}

# Function to calculate average quality
function aveQual(quals, len) {
    sum = 0
    # parse all Q scores
    for (i = 1; i <= len; i++) {
        cmd = "printf '%s' \"" substr(quals, i, 1) "\" | od -A n -t u1"
        cmd | getline q
        close(cmd)
        sum += 10^((q-33) / -10)  # Adjust for ASCII offset of FastQ quality scores
    }
    avg_prob = sum / len
    avg_phred = -10 * log(avg_prob) / log(10)
    return avg_phred
}

# Function to directly subtract 33 from ASCII value
function ord(c) {
    return int(sprintf("%d", -33 + c))
}

# round decimal numbers
function round_decimal(number, decimals) {
    return sprintf("%." decimals "f", number)
}

# Main AWK program
BEGIN {
    FS = "\t"
    print "readid,meanq,length,gc"
}

# Process FastQ records
/^@/ {
    # read 4 rows = 1 record
    header=$1
    getline seq
    getline sep
    getline qual

    # Extract the first field from the space-separated header and remove leading '@'
    split(header, header_array, " ")
    read_id = header_array[1]
    gsub(/^@/, "", read_id)

    # Calculate length of the second row
    seq_length = length(seq)

    # Calculate average quality using aveQual function
    average_quality = aveQual(qual, seq_length)

    # Calculate percentage of G+C letters in the sequence
    gc_count = gsub(/[GCgc]/, "", seq)
    gc_percentage = round_decimal(gc_count / seq_length, 2)

    # Print the results in comma-separated format
    print read_id "\t" average_quality "\t" seq_length "\t" gc_percentage
}