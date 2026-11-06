@loop inner (int n) : k { for (int k = 0; k < n; k++) {...} }
@loop given_too_late (int n) : k { inner (n) { } {...} }
@loop wrong_count (int n) { inner (n, n) { {...} } }
@loop returns (int n) { inner (n) { return; } {...} }
@loop itself (int n) { inner (n) { itself (n) { } } {...} }
@loop hides_k (int n) { inner (n) { {...} } }
void f (void) { hides_k (2) { (void) k; } inner (1) inner (1, 2) { } }
