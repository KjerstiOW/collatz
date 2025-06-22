import math
from fractions import Fraction
from sympy.functions import totient

def check_parity(n: int):
    """
    Throws an error if `n` is odd.

    Parameters:
        n (int): The integer to check.

    Returns:
        bool: True if `n` is odd.
    
    Raises:
        ValueError: If `n` is not odd.
    """

    if n % 2 == 0:
        raise ValueError(f"{n} must be an odd integer")
    return True

def discrete_log(base: int, result: int, modulo: int):
    """
    Computes the discrete logarithm.
    Solves for the smallest non-integer `k` such that:
        base^k ≡ result (mod modulo)

    Parameters:
        base (int): The base of the exponent.
        result(int): The target value of `base^k % modulo`.
    
    Returns:
        (int | None): The smallest non-negative integer `k` satisfying the congruence. If no `k` is found, returns `None`.
    
    Note:
        This uses brute-force search, so it's inefficient for large moduli.
    """

    for k in range(modulo):
        if pow(base, k, modulo) == result:
            return k
        
    return None

def sum_exponents(lst: list[int]):
    """
    Computes the cumulative summations for the encoded Collatz list.
    
    Parameters:
        lst (list[int]): The encoded Collatz list.
    
    Returns:
        list[int]: The cumulative summations (s_1, s_2, ..., s_m).
    
    Example:
        >>> sum_exponents([2, 3, 4])
        [9, 5, 2]
    """
    exponents = [sum(lst[::-1][k:]) for k in range(len(lst))]

    return exponents

def encode_integer(n: int):
    """
    Encode an odd integer into its inverse Collatz list.

    Parameters:
        n (int): An odd positive integer to encode.
    
    Returns:
        list[int]: The inverse Collatz encoding of `n`.
    
    Raises:
        ValueError: If `n` is not odd.
    
    Example:
        >>> encode_integer(7)
        [1, 1, 2, 3, 4]
    """
    check_parity(n)

    lst = []

    while n != 1:
        if n % 2 == 0:
            lst[-1] += 1
            n //= 2
        else:
            lst.append(0)
            n = 3*n + 1
    
    return lst

def leading_term_from_list(lst: list[int]):
    """
    Computes the leading term of the algebraic Collatz expression.

    Parameters:
        lst (list[int]): The Collatz exponent list.
    
    Returns:
        Fraction: The leading term of the algebraic Collatz expression.
    
    Example:
        >>> leading_term_from_list([2, 3, 4])
        Fraction(512, 27)
    """
    exponents = sum_exponents(lst)
    return Fraction(2**exponents[0], 3**len(lst))

def leading_term_from_int(n: int):
    """
    Computes the leading term of the given odd integer.

    Parameters:
        n (int): An odd integer.
    
    Returns:
        Fraction: The leading term of `n`.

    Raises:
        ValueError: If `n` is not odd.
    
    Example:
        >>> leading_term_from_list(7)
        Fraction(512, 27)
    """    
    check_parity(n)

    lst = encode_integer(n)
    return leading_term_from_list(lst)

def remainder_from_list(lst: list[int]):
    """
    Computes the remainder of the encoded Collatz list.

    Parameters:
        lst (list[int]): The encoded Collatz list.
    
    Returns:
        Fraction: The remainder of the encoded Collatz list `lst`.
    
    Example:
        >>> remainder_from_list([2, 3, 4])
        Fraction(65, 27)
    """

    exponents = sum_exponents(lst)
    m = len(exponents)
    R = 0
    for e in range(1, m):
        R += Fraction(2**lst[e], 3**(m-e+1))
    return R + Fraction(1, 3)

def remainder_from_int(n: int):
    """
    Computes the remainder of the given odd integer.

    Parameters:
        n (int): An odd integer.
    
    Returns:
        Fraction: The Collatz remainder of the odd integer `n`.
    
    Raises:
        ValueError: If `n` is not odd.
    
    Example:
        >>> remainder_from_int(7)
       Fraction(65, 27)
    """
    check_parity(n)

    lst = encode_integer(n)
    return remainder_from_list(lst)

def algebraic_collatz(n_values: list[int]):
    """
    Returns the algebraic Collatz expression at the given encoded list.

    Parameters:
        n_values (list[int]): An encoded Collatz list.
    
    Returns:
        Fraction: The numerical value of the encoded Collatz list.
    
    Example:
        >>> algebraic_collatz([2, 3, 4])
        Fraction(7, 1)
    """
    exponents = sum_exponents(n_values)
    m = len(exponents)
    n = Fraction(2**exponents[0], 3**m)

    for e in range(1, m):
        n -= Fraction(2**exponents[e], 3**(m-e+1))
    return n - Fraction(1, 3)

def find_first_n(lst: list[int]):
    """
    Finds `n_1` given `n_2` through `n_m`.

    Parameters:
        lst (list[int]): The encoded Collatz list `[n_m, ..., n_2]`.

    Returns:
        (int | None): The value of `n_1`, None if no value is found.
    
    Example:
        >>> find_first_n([2, 3])
        4
    """
    exponents = sum_exponents(lst)
    m = len(lst)
    modulus = 3**(m+1)
    C = 0

    for i, exp in enumerate(exponents):
        C += 3**(i) * 2**exp % modulus

    C += 3**(m)

    phi = totient(modulus)
    log2C = discrete_log(2, C % modulus, modulus)

    if log2C is None:
        return None

    return (log2C - exponents[0]) % phi

def list_has_extension(lst: list[int]):
    """
    Determines if a list has a Collatz extension.

    Parameters:
        lst (list[int]): The encoded Collatz list.
    
    Returns:
        bool: `True` if `lst` has an extension, `False` if `lst` has no extension.
    
    Example:
        >>> list_has_extension([2, 3, 4])
        True
    """
    if algebraic_collatz(lst) % 3 == 0:
        return False
    return True

def int_has_extension(n: int):
    """
    Determines if an integer has a Collatz extension.

    Parameters:
        n (int): An integer.
    
    Returns:
        bool: `True` if `n` has a list extension, `False` if `n` has no list extension.
    
    Example:
        >>> int_has_extension(7)
        True
    """
    if n % 3 == 0:
        return False
    return True

def get_representative_count(m: int):
    """
    Calculates the representative count of length `m`.

    Parameters:
        m (int): The list length to calculate.
    
    Returns:
        int: The number of representatives of length `m`.
    
    Example:
        >>> get_representative_count(3)
        12
    """

    return 2**(m-1) * 3**((m-1) * (m-2)/2)

def find_m_s1_pairs(n, max_length=None, alpha=1.26):
    """
    Finds all valid (m, s_1) pairs for the given integer.

    Parameters:
        n (int): A positive odd integer to analyze.
        max_length (int, optional): Maximum number of odd steps to check. Defaults to floor(14 * log2(n)).
        alpha(float, optional): The upper bound on alpha. Defaults to 1.26.
    
    Returns:
        list[tuple[int, int]]: A list of all (m, s_1) integer pairs that satisfy the bound.
    
    Raises:
        ValueError: If `n` is not odd.
        
    Example:
        >>> find_m_s1_pairs(7)
        [(2, 6), (5, 11), (7, 14), (12, 22), (14, 25), (17, 30), (19, 33), (22, 38), (24, 41), (29, 49), (31, 52), (34, 57), (36, 60)]
    """
    check_parity(n)

    if max_length == None:
        max_length = math.floor(14 * math.log2(n))

    nums = range(1, max_length+1)

    possible_values = []
    for i in nums:
        lower_bound = n * 3**i
        upper_bound = n * alpha * 3**i
        j = 1
        while True:
            if (lower_bound <= 2**j and 2**j <= upper_bound):
                possible_values.append((i, j))
            if (2**j > upper_bound):
                break
            j += 1
    
    return possible_values

def v2(n: int) -> int:
    """
    Computes the 2-adic valuation.

    Parameters:
        n (int): A nonzero integer.
    
    Returns:
        int: The largest exponent `v` such that `2^v` divides `n`.
    
    Example:
        >>> v2(8)
        3
    """
    if n == 0:
        raise ValueError("v₂(0) is infinite / undefined")
    return (n & -n).bit_length() - 1

def find_sm(n, m, s1):
    """
    Computes s_m given an odd integer `n`, list length `m`, and cumulative summation `s1`.

    Parameters:
        n (int): An odd integer.
        m (int): The list length.
        s1 (int): The cumulative summation.

    Returns:
        int: The sum s_m.
    """
    return v2(2**s1 - (3**m * n) - 3**(m-1))

def next_n(n, nm):
    """
    Apply one Collatz step; computes the new odd integer.

    Parameters:
        n (int): Current odd integer.
        nm (int): The number of divisions to apply.
    
    Returns:
        (int | None): The new odd integer if the result is integral, `None` otherwise.

    Example:
        >>> next_n(5, 3)
        1
    """
    candidate = Fraction(3 * n + 1, 2**nm)
    return int(candidate) if candidate.denominator == 1 else None

def generate_exponent_list(n, m, s1):
    """
    Generates a valid Collatz exponent list from the given odd integer `n`, list length `m`, and cumulative sum `s1`.

    Parameters:
        n (int): The target odd integer.
        m (int): The of total odd steps to apply.
        s1 (int): The total number of divisions by 2.
    
    Returns:
        (list[int] | None): The constructed list `(n_m, ..., n_1)` or `None` if the list is invalid.

    Raises:
        ValueError: If `n` is not odd.

    Example:
        >>> generate_exponent_list(7, 5, 11)
        [1, 1, 2, 3, 4]
    """
    check_parity(n)
    lst = []
    current_n, current_m, current_s1 = n, m, s1

    while current_m > 1:
        s_m = find_sm(current_n, current_m, current_s1)
        current_n = next_n(current_n, s_m)
        if current_n is None:
            return None

        lst.append(s_m)
        current_m -= 1
        current_s1 -= s_m
    
    lst.append(current_s1)
    
    if lst[-1] % 2 == 1:
        return None
    
    if lst[-1] % 3 == 0 and len(lst) != 1:
        return None

    return lst
def get_candidate_lists(n: int, max_length=None, alpha=1.26):
    """
    Generates all valid exponent lists (n_m, ..., n_1) that could converge to `n`.

    Parameters:
        n (int): A positive odd integer.
        max_length (int, optional): Maximum number of odd steps to consider. Defaults to 14log2(n).
        alpha(float, optional): Controls the upper bound for the leading term. Defaults to 1.26.
    
    Returns:
        list[int]: A list of valid exponent lists that could converge to `n`.
    
    Raises:
        ValueError: If `n` is not odd.

    Example:
        >>> get_candidate_lists(7)
        [[1, 1, 2, 3, 4]]
    """
    check_parity(n)
    candidates = []
    pairs = find_m_s1_pairs(n,max_length=max_length, alpha=alpha)
    for (m, s1) in pairs:
        lst = generate_exponent_list(n, m, s1)
        if lst is not None:
            candidates.append(lst)
    return candidates