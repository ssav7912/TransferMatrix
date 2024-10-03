#pragma once
#include <vector>
#include <fstream>
#include <sstream>
#include <format>

#define CLAMP(v, vmin, vmax) ((std::max)((std::min)((v), (vmax)), (vmin)))


struct lut_range {

    float min;
    float max;

    lut_range(float _min = std::numeric_limits<float>::quiet_NaN(),
        float _max = std::numeric_limits<float>::quiet_NaN()) {

        min = _min;
        max = _max;
    }

    bool operator!=(const lut_range& other) const {

        return min != other.min || max != other.max;
    }
};

/// <summary>
/// LUT Sampler from Randrianasandra et al. (2022) Transfer-Matrix implementation.
/// </summary>
/// <typeparam name="LutDataType"></typeparam>
/// <typeparam name="LutDimension"></typeparam>
template<int LutDimension, typename LutDataType>
class lut {
private:

    // Data size
    int _size[LutDimension];
    // Data range
    lut_range _range[LutDimension];
    // Data
    std::vector<LutDataType> _data;

public:

    lut() { }
    lut(const std::string& filename) {

        std::ifstream in(filename, std::ios_base::in | std::ios_base::binary);

        // Warn if not loaded
        if (in.bad() || in.fail())
            fprintf(stderr, "Unable to load %dD LUT data file: %s\n", LutDimension, filename.c_str());
        else
            fprintf(stdout, "Loading %dD LUT data file: %s\n", LutDimension, filename.c_str());

        this->load(in);
    }
    lut(const int size[], const lut_range range[]) {

        for (int i = 0; i < LutDimension; ++i) {

            _size[i] = size[i];
            _range[i] = range[i];
        }

        _data.assign(linear_size(), LutDataType());
    }

    const int* size() const { return _size; }
    constexpr int dimensions() const { return LutDimension; }
    const lut_range* range() const { return _range; }

    LutDataType* data() { return _data.data(); }
    const LutDataType* data() const { return _data.data(); }

    std::string size_string() const {

        std::ostringstream oss;

        for (int i = 0; i < LutDimension; ++i) {
            oss << _size[i];
            if (i < (LutDimension - 1))
                oss << "x";
        }

        return oss.str();
    }
    std::string range_string() const {

        std::ostringstream oss;
        for (int i = 0; i < LutDimension; ++i) {

            oss << "[" << _range[i].min << "," << _range[i].max << "]";
            if (i < (LutDimension - 1))
                oss << "x";
        }
        return oss.str();
    }

    int linear_size() const {

        int lsize = 1;
        for (int i = 0; i < LutDimension; ++i)
            lsize *= _size[i];
        return lsize;
    }
    int linear_index_from_grid_index(const int grid_index[]) const {

        int index = grid_index[LutDimension - 1];
        for (int i = LutDimension - 2; i >= 0; --i)
            index = index * _size[i] + grid_index[i];
        return index;
    }

    void load(std::ifstream& in) {

        // Read size
        for (int i = 0; i < LutDimension; ++i)
            in.read((char*)&_size[i], sizeof(int));

        fprintf(stdout, "Loading %dD LUT of dimensions %s\n", LutDimension, size_string().c_str());

        // Read data range (min / max)
        for (int i = 0; i < LutDimension; ++i) {

            in.read((char*)&_range[i].min, sizeof(float));
            in.read((char*)&_range[i].max, sizeof(float));
        }

        fprintf(stdout, "Loading %dD LUT data range: %s\n", LutDimension, range_string().c_str());

        // Read data
        _data.assign(linear_size(), LutDataType());
        in.read((char*)_data.data(), linear_size() * sizeof(LutDataType));
    }

    template<typename... T>
    inline LutDataType range_get(T... vrange_coords) const {

        float	x[LutDimension] = { vrange_coords... };
        int		i[LutDimension];

        for (int k = 0; k < LutDimension; ++k) {

            // Integer index (ensure the index stays in the limits)
            i[k] = (int)floor((_size[k] - 1) * (x[k] - _range[k].min) / (_range[k].max - _range[k].min));
            i[k] = CLAMP(i[k], 0, _size[k] - 1);
        }

        return _data[linear_index_from_grid_index(i)];
    }
    template<typename... T>
    inline LutDataType range_get_interpolate(T... vrange_coords) const {

        float	x[LutDimension] = { vrange_coords... };
        float	a[LutDimension];
        int		i[LutDimension];
        float	alphas[LutDimension];

        for (int k = 0; k < LutDimension; ++k) {

            // Floating point index
            a[k] = (_size[k] - 1) * (x[k] - _range[k].min) / (_range[k].max - _range[k].min);

            // Integer index
            i[k] = (int)floor(a[k]);

            // Ensure the indexes stays in the limits
            i[k] = CLAMP(i[k], 0, _size[k] - 1);


            // Clamp the interpolation weights
            alphas[k] = CLAMP(a[k] - (float)i[k], 0.f, 1.f);
        }

        int index = linear_index_from_grid_index(i);

        // Lookup result
        LutDataType v = LutDataType();

        // For every possible combinaison of index shift per dimension,
        // fetch the value in memory and do linear interpolation.
        // We fetch using shift of 0 and 1.
        //
        //     v(i+di, j+di, k+dk, l+dl),  where dk in [0,1]
        //
        const unsigned int D = (unsigned int)pow(2, LutDimension);
        for (unsigned int d = 0; d < D; ++d) {

            float alpha = 1.0; // Global alpha
            int   cid_s = 0;   // Id shift

            // Evaluate the weight of the sample d which correspond to
            // one for the shifted configuration:
            // The weight is the product of the weights per dimension.
            for (int k = (LutDimension - 1); k >= 0; --k) {

                bool  bitset = ((1 << k) & d);
                float calpha = (bitset) ? alphas[k] : 1.f - alphas[k];

                // Correct the shift to none if we go out of the grid
                if (i[k] + 1 >= _size[k])
                    bitset = false;

                alpha *= calpha;
                cid_s = cid_s * _size[k] + ((bitset) ? 1 : 0);
            }

            v += alpha * _data[index + cid_s];
        }

        return v;
    }
};

typedef lut<2, float> lut2;
typedef lut<3, float> lut3;
typedef lut<4, float> lut4;