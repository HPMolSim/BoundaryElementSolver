# this script is for the energy and potential of a bi-sphere system, based on the matlab code: https://github.com/zcgan/two_sph_energy_by_img_reflection

export BiSphereSystem, bisphere_energy, bisphere_phi

struct BiSphereSystem
    q_1::Float64
    q_2::Float64
    r_1::Float64
    r_2::Float64
    z_1::Float64
    z_2::Float64
    eps_0::Float64
    eps_1::Float64
    eps_2::Float64
    images::Array{Float64, 2}
end

all_fields(bisys::BiSphereSystem) = (bisys.q_1, bisys.q_2, bisys.r_1, bisys.r_2, bisys.z_1, bisys.z_2, bisys.eps_0, bisys.eps_1, bisys.eps_2, bisys.images)

function BiSphereSystem(q_1::Float64, q_2::Float64, r_1::Float64, r_2::Float64, z_1::Float64, z_2::Float64, eps_0::Float64, eps_1::Float64, eps_2::Float64)
    images = bisphere_images(q_1, q_2, r_1, r_2, z_1, z_2, eps_0, eps_1, eps_2)
    return BiSphereSystem(q_1, q_2, r_1, r_2, z_1, z_2, eps_0, eps_1, eps_2, images)
end

function jacobi_ab(n::Int, alpha::Float64, beta::Float64, a::Float64, b::Float64)

    x, w = jacobi(n, alpha, beta)
    x = (b-a)/2 * x .+ (a+b)/2
    w = ((b-a)/2)^(1 + alpha + beta) * w
    return x, w
end

function bisphere_energy(bisys::BiSphereSystem)

    q_1, q_2, r_1, r_2, z_1, z_2, eps_0, eps_1, eps_2, img = all_fields(bisys)

    energy = 0.0
    # Energy calculation
    dx = 0.0 - 0.0
    dy = 0.0 - 0.0
    dz = z_1 - z_2
    dr = sqrt(dx^2 + dy^2 + dz^2)
    # energy += q_1 * q_2 / dr / (4π * eps_0)

    src = zeros(Float64, 2, 4)  # 2 charges, each has its x,y,z coordinates and charge q
    src[1, :] = [0.0, 0.0, z_1, q_1]
    src[2, :] = [0.0, 0.0, z_2, q_2]

    # Source-image interaction (i.e., the polarization part)
    # For first sphere
    for j in 1:size(img, 1)
        dx = img[j, 1] - src[1, 1]
        dy = img[j, 2] - src[1, 2]
        dz = img[j, 3] - src[1, 3]
        dr = sqrt(dx^2 + dy^2 + dz^2)
        if dr > r_1
            energy += 0.5 * q_1 * img[j, 4] / dr / (4π * eps_0)
        end
    end

    # For second sphere
    for j in 1:size(img, 1)
        dx = img[j, 1] - src[2, 1]
        dy = img[j, 2] - src[2, 2]
        dz = img[j, 3] - src[2, 3]
        dr = sqrt(dx^2 + dy^2 + dz^2)
        if dr > r_2
            energy += 0.5 * q_2 * img[j, 4] / dr / (4π * eps_0)
        end
    end

    return energy
end

function bisphere_phi(x, y, z, bisys::BiSphereSystem)

    q_1, q_2, r_1, r_2, z_1, z_2, eps_0, eps_1, eps_2, images = all_fields(bisys)

    dr_1 = sqrt(x^2 + y^2 + (z - z_1)^2)
    dr_2 = sqrt(x^2 + y^2 + (z - z_2)^2)
    # @assert dr_1 ≥ r_1 && dr_2 ≥ r_2 "Position is inside the sphere"

    phi = q_1 / dr_1 + q_2 / dr_2

    for j in 1:size(images, 1)
        dx = x - images[j, 1]
        dy = y - images[j, 2]
        dz = z - images[j, 3]
        dr = sqrt(dx^2 + dy^2 + dz^2)

        # check if the point is outside the spheres
        # dr_1 = sqrt(x^2 + y^2 + (z - z_1)^2)
        # dr_2 = sqrt(x^2 + y^2 + (z - z_2)^2)

        # if ((dr ≥ r_1) || (dr ≈ r_1)) && ((dr ≥ r_2) || (dr ≈ r_2))
        #     phi += images[j, 4] / dr
        # end

        phi += images[j, 4] / dr
    end

    return phi / (4π * eps_0)
end

function bisphere_images(q_1, q_2, r_1, r_2, z_1, z_2, eps_0, eps_1, eps_2)
    # Set the source and image charge array
    src = zeros(Float64, 2, 4)  # 2 charges, each has its x,y,z coordinates and charge q
    img_max = 1000000  # the maximum number of images allowed
    img = zeros(Float64, img_max, 4)  # the array storing the image location and charges
    img_counter = 1  # integer for counting the number of images

    # Initialize src array
    src[1, :] = [0.0, 0.0, z_1, q_1]
    src[2, :] = [0.0, 0.0, z_2, q_2]

    # Initialize the image generation process
    nlevel = 6  # image series is reflected nlevel times
    order = zeros(Int, nlevel)  # order is the number of discretized Gauss-Jacobi quadrature points
    order[1] = 5  # first level
    order[2] = 4  # second level
    for i in 3:nlevel
        order[i] = 3  # higher levels
    end

    # Generate the first level images due to the two point sources
    for i in 1:2
        pos = src[i, :]
        img_counter = generate_img!(pos, order[1], r_1, r_2, z_1, z_2, eps_0, eps_1, eps_2, img, img_counter)
    end

    # Generate the higher level images due to reflection
    img_ind_start = 1
    img_ind_end = img_counter - 1
    for i in 2:nlevel
        for j in img_ind_start:img_ind_end
            pos = img[j, :]
            img_counter = generate_img!(pos, order[i], r_1, r_2, z_1, z_2, eps_0, eps_1, eps_2, img, img_counter)
        end
        img_ind_start = img_ind_end + 1
        img_ind_end = img_counter - 1
    end

    return img[1:img_counter - 1, :]
end

function generate_img!(pos, order, r_1, r_2, z_1, z_2, eps_0, eps_1, eps_2, img, img_counter)
    # Constants related to eps_i
    eps_i = [eps_1, eps_2]
    lamda = eps_0 ./ (eps_0 .+ eps_i)
    beta = lamda .- eps_0
    gamma = (eps_i .- eps_0) ./ (eps_0 .+ eps_i)

    # Image charge in the first sphere
    dx = pos[1] - 0.0
    dy = pos[2] - 0.0
    dz = pos[3] - z_1
    dr = sqrt(dx^2 + dy^2 + dz^2)

    if dr > r_1
        rk = r_1^2 / dr
        rim = rk / dr
        img[img_counter, 1] = rim * dx + 0.0
        img[img_counter, 2] = rim * dy + 0.0
        img[img_counter, 3] = rim * dz + z_1
        img[img_counter, 4] = -gamma[1] * r_1 * pos[4] / dr
        img_counter += 1

        # Get Gauss quadrature weights and locations
        alpha = 0.0  # always zero
        a = 0.0
        b = rk  # a is always 0
        x, w = jacobi_ab(order, alpha, beta[1], a, b)
        for i in 1:order
            img[img_counter, 1] = 0 + x[i] * dx / dr
            img[img_counter, 2] = 0 + x[i] * dy / dr
            img[img_counter, 3] = z_1 + x[i] * dz / dr
            img[img_counter, 4] = w[i] * gamma[1] * lamda[1] * rk^(1 - lamda[1]) * pos[4] / r_1
            img_counter += 1
        end
    end

    # Image charge in the second sphere
    dx = pos[1] - 0.0
    dy = pos[2] - 0.0
    dz = pos[3] - z_2
    dr = sqrt(dx^2 + dy^2 + dz^2)

    if dr > r_2
        rk = r_2^2 / dr
        rim = rk / dr
        img[img_counter, 1] = rim * dx + 0.0
        img[img_counter, 2] = rim * dy + 0.0
        img[img_counter, 3] = rim * dz + z_2
        img[img_counter, 4] = -gamma[2] * r_2 * pos[4] / dr
        img_counter += 1

        # Get Gauss quadrature weights and locations
        alpha = 0.0  # always zero
        a = 0.0
        b = rk  # a is always 0
        x, w = jacobi_ab(order, alpha, beta[2], a, b)
        for i in 1:order
            img[img_counter, 1] = 0 + x[i] * dx / dr
            img[img_counter, 2] = 0 + x[i] * dy / dr
            img[img_counter, 3] = z_2 + x[i] * dz / dr
            img[img_counter, 4] = w[i] * gamma[2] * lamda[2] * rk^(1 - lamda[2]) * pos[4] / r_2
            img_counter += 1
        end
    end

    return img_counter
end