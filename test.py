import math
from decimal import Decimal, getcontext
from mpmath import radians, cos, sin

# Set the precision for Decimal operations (e.g., 50 decimal places)
getcontext().prec = 50

# Earth's radius in kilometers
R = Decimal(6371.0)

# Function to convert degrees to radians using mpmath for higher precision
def deg_to_rad(degrees):
    return radians(degrees)

# Function to convert spherical coordinates (Longitude, Latitude, Depth) to Cartesian (X, Y, Z)
def spherical_to_cartesian(lon, lat, depth):
    # Convert longitude and latitude from degrees to radians using mpmath
    lon_rad = deg_to_rad(lon)
    lat_rad = deg_to_rad(lat)
    
    # Calculate Cartesian coordinates (X, Y, Z) with higher precision using mpmath
    X = R * cos(lat_rad) * cos(lon_rad)
    Y = R * cos(lat_rad) * sin(lon_rad)
    Z = R * sin(lat_rad) - Decimal(depth)  # Adjust for depth

    return X, Y, Z

# Function to read an XYZ file, convert its spherical coordinates to Cartesian, and return the result
def convert_xyz_file(input_file, output_file):
    with open(input_file, 'r') as infile:
        with open(output_file, 'w') as outfile:
            for line in infile:
                # Read each line from the XYZ file
                lon, lat, depth = map(float, line.split())
                
                # Convert the spherical coordinates to Cartesian with high precision
                X, Y, Z = spherical_to_cartesian(lon, lat, depth)
                
                # Write the Cartesian coordinates to the output file with high precision
                outfile.write(f"{X:.15f} {Y:.15f} {Z:.15f}\n")
    
    print(f"Conversion complete. Results saved to {output_file}")

# Example Usage
input_file = 'selatsunda.xyz'  # Replace with the path to your input XYZ file
output_file = 'selatsundac2.xyz'  # The path where the Cartesian coordinates will be saved
convert_xyz_file(input_file, output_file)