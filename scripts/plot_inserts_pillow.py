from PIL import Image, ImageDraw

def text_to_bitmap(text_file, output_image):
    colors = {
        'A': 'green',
        'C': 'red',
        'G': 'blue',
        'T': 'yellow',
        'a': 'green',
        'c': 'red',
        'g': 'blue',
        't': 'yellow',
        '-': 'black',
        '>': 'black'
    }
    
    # Read lines from the text file
    with open(text_file, 'r') as file:
        lines = file.readlines()

    lines = [line for line in lines if len(line) <= 300]

    # Determine image dimensions
    width = max(len(line) for line in lines)*12
    height = len(lines)

    width = 300*12

    # Create a new image with a white background
    image = Image.new('RGB', (width, height), 'white')

    # Draw each pixel based on the text content
    draw = ImageDraw.Draw(image)
    for y, line in enumerate(lines):
        for x, char in enumerate(line):
            if char != ' ' and char != '\n':
                for offset in range(12):
                    draw.point((x*12+offset, y), fill=colors[char])
        draw.point((15*12, y), fill='grey')
        draw.point((48*12, y), fill='grey')
        draw.point((111*12, y), fill='grey')
        draw.point((249*12, y), fill='grey')

    # Save the image
    image.save(output_image)

if __name__ == "__main__":
    # text_file_path = "../notebooks/fwd_ins_ext_00_05.txt"
    # output_image_path = "../notebooks/fwd_in_ext_00_05_bitmap.png"

    # text_to_bitmap(text_file_path, output_image_path)

    # text_file_path = "../notebooks/rev_ins_ext_00_05.txt"
    # output_image_path = "../notebooks/rev_in_ext_00_05_bitmap.png"

    # text_to_bitmap(text_file_path, output_image_path)

    text_file_path = "../notebooks/fwd_ins_ext_00_67.txt"
    output_image_path = "../notebooks/fwd_in_ext_00_67_bitmap.png"

    text_to_bitmap(text_file_path, output_image_path)

    text_file_path = "../notebooks/rev_ins_ext_00_67.txt"
    output_image_path = "../notebooks/rev_in_ext_00_67_bitmap.png"

    text_to_bitmap(text_file_path, output_image_path)
