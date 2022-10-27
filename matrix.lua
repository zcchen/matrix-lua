local matrix ={}

local mat_obj = {
    data = {}
}
local mat_meta = {}

local function deep_copy(object)
    local lookup_table = {}
    local function _copy(obj)
        if type(obj) ~= "table" then
            return obj
        elseif lookup_table[obj] then
            return lookup_table[obj]
        end
        local new_table = {}
        lookup_table[obj] = new_table
        for key, value in pairs(obj) do
            new_table[_copy(key)] = _copy(value)
        end
        return setmetatable(new_table, getmetatable(obj))
    end
    return _copy(object)
end

---------------- init the matrix --------------------
local function mat_init()
    local mat = {}
    return deep_copy(mat_obj)
end

function matrix.zeros(row, column)
    local mat = mat_init()
    for i = 1, row do
        mat.data[i] = {}
        for j = 1, column do
            mat.data[i][j] = 0
        end
    end
    setmetatable(mat, mat_meta)
    return mat
end

function matrix.I(m)
    local mat = matrix.zeros(m, m)
    for i = 1, m do
        mat.data[i][i] = 1
    end
    return mat
end

function matrix.new(row, ...)
    local rows = {row, ...}
    local column = #row
    local ret = matrix.zeros(#rows, column)
    for i = 1, #rows, 1 do
        local r = rows[i]
        if column ~= #r then
            print("[ERROR] At least one row has different column than others.")
            return nil
        end
        for j = 1, column, 1 do
            local c = r[j]
            if type(c) ~= "number" then
                print("[ERROR] At least one element is NOT type number.")
                return nil
            else
                ret.data[i][j] = c
            end
        end
    end
    return ret
end
-----------------------------------------------------

------------- base function of the matrix -----------
function mat_obj:size()
    local row, column = #self.data, #self.data[1]
    return row, column
end

function mat_obj:print()
    local row, column = self:size()
    for i = 1, row do
        for j =1, column do
            io.write(self.data[i][j])
            io.write(" ")
        end
        print("")
    end
end

function mat_meta:__add(mat)
    if getmetatable(mat) ~= mat_meta then
        print("[ERROR] matrix can add / sub matrix ONLY.")
        return nil
    end
    local self_row, self_column = self:size()
    local mat_row, mat_column = mat:size()
    if self_row ~= mat_row or self_column ~= mat_column  then
        print(string.format(
                "[ERROR] matrix <%dx%d> cannot add / sub with <%dx%d>",
                self_row, self_column, mat_row, mat_column
            )
        )
        return nil
    end
    local ret = matrix.zeros(self_row, self_column)
    for i = 1, self_row, 1 do
        for j = 1, self_column, 1 do
            ret.data[i][j] = self.data[i][j] + mat.data[i][j]
        end
    end
    return ret
end

function mat_meta:__mul(mat_or_num)
    if getmetatable(mat_or_num) == mat_meta then
        local self_row, self_column = self:size()
        local mat_row, mat_column = mat_or_num:size()
        if self_column ~= mat_row  then
            print(string.format(
                    "[ERROR] matrix <%dx%d> cannot multiply with <%dx%d>",
                    self_row, self_column, mat_row, mat_column
                )
            )
            return nil
        end
        local ret = matrix.zeros(self_row, mat_column)
        for i = 1, self_row, 1 do
            for j = 1, mat_column, 1 do
                for k = 1, self_column, 1 do
                    ret.data[i][j] = ret.data[i][j] + self.data[i][k] * mat_or_num.data[k][j]
                end
            end
        end
        return ret
    elseif type(mat_or_num) == "number" then
        local ret = deep_copy(self)
        local self_row, self_column = self:size()
        for i = 1, self_row, 1 do
            for j = 1, self_column, 1 do
                ret.data[i][j] = ret.data[i][j] * mat_or_num
            end
        end
        return ret
    else
        print("[ERROR] matrix can add matrix or number ONLY.")
        return nil
    end
end

function mat_meta:__sub(mat)
    return self + mat * -1
end

function mat_meta:__eq(mat)
    if getmetatable(mat) ~= mat_meta then
        print("[ERROR] matrix can compare with matrix ONLY.")
        return nil
    end
    local self_row, self_column = self:size()
    local mat_row, mat_column = mat:size()
    if self_row ~= mat_row or self_column ~= mat_column  then
        return false
    else
        for i= 1, self_row, 1 do
            for j = 1, self_column, 1 do
                if self.data[i][j] ~= mat.data[i][j] then
                    return false
                end
            end
        end
        return true
    end
end

function mat_obj:T()
    local self_row, self_column = self:size()
    local ret = matrix.zeros(self_column, self_row)
    for i = 1, self_column, 1 do
        for j = 1, self_row, 1 do
            ret.data[i][j] = self.data[j][i]
        end
    end
    return ret
end

-- ref: https://rosettacode.org/wiki/LU_decomposition#JavaScript
function mat_obj:pivot()
    local self_row, self_column = self:size()
    if self_row ~= self_column then
        print(string.format(
            "[ERROR] matrix <%dx%d> does NOT have pivotized matrix",
            self_row, self_column)
        )
        return nil
    end
    local ret = matrix.I(self_row)
    for i = 1, self_column, 1 do
        local maxm = self.data[i][i]
        local curr_row = i
        for j = i, self_row, 1 do
            if self.data[j][i] > maxm then
                maxm = m[j][i]
                row = j
            end
        end
        if i ~= curr_row then   -- switch the lines
            local tmp = ret.data[i]
            ret.data[i] = ret.data[row]
            ret.data[row] = tmp
        end
    end
    return ret
end

function mat_obj:LUdecompose()
    -- PLU decompostion solution
    -- ref: https://en.wikipedia.org/wiki/LU_decomposition
    local self_row, self_column = self:size()
    if self_row ~= self_column then
        print(string.format(
            "[ERROR] matrix <%dx%d> cannot be L-U decomposed.",
            self_row, self_column)
        )
        return nil
    end
    local L = matrix.I(self_row)
    local U = matrix.zeros(self_row, self_column)
    for i = 1, self_row, 1 do
        for k = i, self_column, 1 do
            local sum = 0
            for j = 1, i, 1 do
                sum = sum + (L.data[i][j] * U.data[j][k])
            end
            U.data[i][k] = self.data[i][k] - sum
        end
        for k = i, self_column, 1 do
            if i == k then
                L.data[i][i] = 1;
            else
                local sum = 0
                for j = 1, i, 1 do
                    sum = sum + L.data[k][j] * U.data[j][i]
                end
                L.data[k][i] = (self.data[k][i] - sum) / U.data[i][i]
            end
        end
    end
    return L, U
end

-- since A = L * U, thus: A^-1 = U^-1 * L^-1
-- matrix U^-1 & L^-1 can be caculated by Elementary row operations
function mat_obj:inverse()
    -- ref: https://stackoverflow.com/questions/420612/is-there-around-a-straightforward-way-to-invert-a-triangular-upper-or-lower-ma
    local esp = 0.0001
    local L, U = self:LUdecompose()
    local n, _ = self:size()
    for i = 1, n, 1 do
        if math.abs(L.data[i][i]) <= esp or math.abs(U.data[i][i]) <= esp then
            print("[ERROR] this matrix cannot be inversed().")
            return nil
        end
    end
    local l_inverse = matrix.zeros(n, n)
    for i = 1, n, 1 do
        l_inverse.data[i][i] = 1 / L.data[i][i]
        for j = 1, i - 1, 1 do
            local s = 0
            for k = j, i - 1, 1 do
                s = s + L.data[i][k] * l_inverse.data[k][j]
            end
            l_inverse.data[i][j] = -s * l_inverse.data[i][i]
        end
    end
    local U_T = U:T()       -- U^-1 = (U^T) ^ -1
    local ut_inverse = matrix.zeros(n, n)
    for i = 1, n, 1 do
        ut_inverse.data[i][i] = 1 / U_T.data[i][i]
        for j = 1, i - 1, 1 do
            local s = 0
            for k = j, i - 1, 1 do
                s = s + U_T.data[i][k] * ut_inverse.data[k][j]
            end
            ut_inverse.data[i][j] = -s * ut_inverse.data[i][i]
        end
    end
    local u_inverse = ut_inverse:T()
    return u_inverse * l_inverse
end
-----------------------------------------------------

return matrix
