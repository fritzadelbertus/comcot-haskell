import TypeModule

-- allGrid:: Layer -> Layer


jnz:: Layer -> Layer -> Layer
jnz lo la
    | sc_option la == 0 = oldZCoupling lo la
    | sc_option la == 1 = newZCoupling lo la

-- oldZCoupling:: Layer -> Layer -> Layer
-- oldZCoupling
--     | 