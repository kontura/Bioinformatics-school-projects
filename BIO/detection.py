import cv2
import numpy 
import sys
import math
from skimage import measure

def display_images(img1, img2):
    img1 = resize(img1)
    img2 = resize(img2)
    cv2.imshow("in", img1)
    cv2.moveWindow("in", 200, 290)
    cv2.imshow('original', img2)
    cv2.moveWindow("original", 800, 290)
    cv2.waitKey(0)
    cv2.destroyAllWindows()
    cv2.imwrite("img2.jpg", img2)
    return

def resize(img, size=500):
    width, height = img.shape[:2]
    return cv2.resize(img, (size, int(math.floor((width/height)*size))))

def load_mask():
    if len(sys.argv) > 2:
        mask = cv2.imread(sys.argv[2])
        if mask is None:
            print ("mask not loaded")
            quit()
        else:
            return mask
    else:
        return None

def load_image():
    if len(sys.argv) > 1:
        img = cv2.imread(sys.argv[1])
        if img is None:
            print ("img not loaded")
            quit()
        else:
            return img
    else:
        print ("Usage: detection.py <image_to_proccess> <image_mask>")
        quit()


def read_images():
    mask = load_mask()
    img = load_image()
    return img, mask

def use_mask(mask, dst):
    mask = cv2.cvtColor(mask, cv2.COLOR_BGR2GRAY)
    dst = cv2.bitwise_and(dst, dst, mask = mask)
    return dst

def crop_image(img):
    dst = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    _, thresh = cv2.threshold(dst, 1, 255, cv2.THRESH_BINARY)
    im2, contours, hierarchy = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    x,y,w,h = cv2.boundingRect(im2)
    crop = img[y:y+h, x:x+w]
    return crop

def get_green_channel(img):
    b,g,r = cv2.split(img)
    #g = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    #as of now i am returing GRAYSCALE !! TODO
    return g

def subtract_background_approximation(img):
    bg_approx = cv2.medianBlur(img, 105)
    dst2 = cv2.absdiff(img , bg_approx)
    return dst2

def shade_correction(image):
    dst = image
    dst = subtract_background_approximation(dst)
    #dst = cv2.normalize(image, dst, alpha=0, beta=1, norm_type=cv2.NORM_MINMAX, dtype=cv2.CV_32F)
    return dst

def reduce_noise(image):
    dst = cv2.erode(image, None, iterations=4)
    dst = cv2.dilate(dst, None, iterations=4)
    return dst

def build_verical_kernel():
    kernel = numpy.full((60,60), -1)
    x = numpy.full((60), -2)
    kernel[:,25]=x
    kernel[:,26]=x
    kernel[:,27]=x
    kernel[:,28]=x
    kernel[:,29]=x

    kernel[:,30]=x
    kernel[:,31]=x
    kernel[:,32]=x
    kernel[:,33]=x
    kernel[:,34]=x
    kernel = kernel.astype(float)
    return kernel

def assign_diagonal_from(x_start, y_start, matrix, value):
    width, height = matrix.shape
    y = y_start
    for x in range (x_start, width):
        if (y < height):
            matrix[x,y] = value
            y += 1
        else:
            return matrix
    return matrix


def build_diagonal_kernel():
    kernel = numpy.full((60,60), -9)
    value = 4
    kernel = assign_diagonal_from(0, 0, kernel, value)
    kernel = assign_diagonal_from(1, 0, kernel, value)
    kernel = assign_diagonal_from(0, 1, kernel, value)
    kernel = assign_diagonal_from(2, 0, kernel, value)
    kernel = assign_diagonal_from(0, 2, kernel, value)
    return kernel


def build_kernels():
    kernels = []
    kernels.append(build_verical_kernel())
    kernels.append(numpy.rot90(kernels[0]))
    kernels.append(build_diagonal_kernel())
    kernels.append(numpy.rot90(kernels[2]))
    numpy.set_printoptions(linewidth=200)
    #print(kernels[2])
    #print(numpy.sum(kernels[2]))
    return kernels

def delete_lines_with_hought(image):
    minLineLength = 1
    maxLineGap = 50
    lines = cv2.HoughLinesP(image,1,numpy.pi/180,100,minLineLength,maxLineGap)
    for line in lines:
        for x1,y1,x2,y2 in line:
            cv2.line(img,(x1,y1),(x2,y2),(0,255,0),10)
    return image

def convolution_with_kernels(image):
    kernels = build_kernels()
    applied_kernels = []
    applied_kernels.append(cv2.filter2D(image, cv2.CV_64F, kernels[3]))
    display_images(applied_kernels[0], image)
    return image

def explore_line(x_to, y_to, output, visited, image):
    output[x_to][y_to] = 1
    visited[x_to][y_to] = 1
    new_x = x_to+1
    new_y = y_to+1
    if ((visited[new_x][new_y] != 1) and (image[new_x][new_y] != 0)):
        explore_line(new_x, new_y, output, visited, image)
    new_x = x_to+0
    new_y = y_to+1
    if ((visited[new_x][new_y] != 1) and (image[new_x][new_y] != 0)):
        explore_line(new_x, new_y, output, visited, image)
    new_x = x_to+1
    new_y = y_to+0
    if ((visited[new_x][new_y] != 1) and (image[new_x][new_y] != 0)):
       explore_line(new_x, new_y, output, visited, image)
    new_x = x_to-1
    new_y = y_to-1
    if ((visited[new_x][new_y] != 1) and (image[new_x][new_y] != 0)):
        explore_line(new_x, new_y, output, visited, image)
    new_x = x_to-1
    new_y = y_to+0
    if ((visited[new_x][new_y] != 1) and (image[new_x][new_y] != 0)):
        explore_line(new_x, new_y, output, visited, image)
    new_x = x_to+0
    new_y = y_to-1
    if ((visited[new_x][new_y] != 1) and (image[new_x][new_y] != 0)):
        explore_line(new_x, new_y, output, visited, image)
    new_x = x_to-1
    new_y = y_to+1
    if ((visited[new_x][new_y] != 1) and (image[new_x][new_y] != 0)):
        explore_line(new_x, new_y, output, visited, image)
    new_x = x_to+1
    new_y = y_to-1
    if ((visited[new_x][new_y] != 1) and (image[new_x][new_y] != 0)):
        explore_line(new_x, new_y, output, visited, image)

def custom_line_searching_by_count(image):
    image = reduce_noise(image)
    smaller_image = resize(image, 300)
    width, height = smaller_image.shape[:2]
    visited = numpy.full((width, height), 0)
    output = numpy.full((width, height), 0)
    for x in range (0, width):
        for y in range (0, height):
            if ((smaller_image[x][y] != 0) and (visited[x][y] != 1)):
                output_to_be_merged = numpy.full((width, height), 0)
                explore_line(x, y, output_to_be_merged, visited, image)
                if (numpy.sum(output_to_be_merged) > 5):
                    display_images(output_to_be_merged, visited)


    display_images(smaller_image, visited)
    return image

def line_extraction_by_border_lenght(image):
    labels = measure.label(image, neighbors=8, background=0)
    mask = numpy.zeros(image.shape, dtype="uint8")
    for label in numpy.unique(labels):
        if label == 0:
            continue
        label_mask = numpy.zeros(image.shape, dtype="uint8")
        label_mask[labels == label] = 255
        cont_detect = numpy.copy(label_mask)
        im2, countours, hierarchy = cv2.findContours(cont_detect, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

        #print(numpy.sum(label_mask)/len(countours[0]))
        #display_images(mask, label_mask)
        if (len(countours[0]) > 300):
            mask = cv2.add(mask, label_mask)
    return mask


def line_extraction(image):
    out = line_extraction_by_border_lenght(image)

    #image = delete_lines_with_hought(image)
    #image = convolution_with_kernels(image)

    #pdb.set_trace()
    #out = custom_line_searching_by_count(image)
    return out

def apply_preprocessing(img, mask):
    if mask is not None:
        img = use_mask(mask, img)
    dst = crop_image(img)
    dst = resize(dst, 1479)
    resized = dst


    #working with green channel, givers better contrast
    dst = get_green_channel(dst)

    #redukce sumu
    #dst = cv2.GaussianBlur(dst, (3,3), 0)
    dst = cv2.medianBlur(dst, 3)

    adaptive_histogram_equalization = cv2.createCLAHE(clipLimit=3.0, tileGridSize=(2,2))
    dst = adaptive_histogram_equalization.apply(dst)

    dst = shade_correction(dst)

    return dst, resized

def vessel_segmentation(dst):
    _, dst = cv2.threshold(dst, 45, 255, cv2.THRESH_BINARY)
    #dst = cv2.adaptiveThreshold(dst, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY, 51, 4)
    lines = line_extraction(cv2.dilate(dst, None, iterations=9))
    #lines = cv2.dilate(lines, None, iterations=19)
    display_images(dst, lines)
    dst = numpy.subtract(dst, lines)

    #dst = reduce_noise(dst)
    return dst

def print_image_as_text(img):
    width, height = img.shape[:2]
    for x in range (0, width):
        for y in range (0, height):
            print(x,y, ":", img[x,y])


img, mask = read_images()
output, resized = apply_preprocessing(img, mask)
output = vessel_segmentation(output)
#print_image_as_text(output)
display_images(resized, output)
