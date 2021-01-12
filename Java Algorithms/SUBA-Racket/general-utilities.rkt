#lang racket
; General utilities
; Dana Angluin
; 9/27/20

(provide entry entry-key entry-value lookup-all lookup
         member? set-diff set-union set-intersect fast-set-intersect set-subset?
         c-product all-subsets all-assigns
         sort-symbols
         for-all? exists?
         pick pick-k random-subset random-sorted-list
         select-index-probs random-select-counts)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; utilities

; A table is a list of entries: each entry has a key and value

(struct entry (key value) #:transparent)

; Look up all the elements of a table that match the given key and
; return a list of all the resulting values with duplicates removed.

(define (lookup-all key table)
  (remove-duplicates
   (map entry-value
        (filter (lambda (item)
                  (equal? key (entry-key item)))
                table))))

; lookup, returning first value that matches key, or #f if none

(define (lookup key table)
  (let ((result (lookup-all key table)))
    (if (empty? result)
        #f
        (first result))))


; membership predicate
(define (member? item lst)
  (if (member item lst) #t #f))

; set difference
(define (set-diff lst1 lst2)
  (filter (lambda (item) (not (member item lst2)))
          lst1))

; set union (with remove-duplicates)
(define (set-union lst1 lst2)
  (remove-duplicates (append lst1 lst2)))

; set intersection
(define (set-intersect lst1 lst2)
  (filter (lambda (item) (member item lst2))
          lst1))

; faster set intersect
(define (fast-set-intersect lst1 lst2)
;  (displayln "fast-set-intersect")
  (let* [(ht (make-hash))]
    (for-each
     (lambda (item)
       (hash-set! ht item #t))
     lst2)
    (filter
     (lambda (item)
       (hash-ref ht item #f))
     lst1)))

; subset
(define (set-subset? lst1 lst2)
  (for-all? (lambda (item) (member item lst2))
          lst1))

; Cartesian product
(define (c-product lst1 lst2)
  (apply append
         (map (lambda (item1)
                (map (lambda (item2)
                       (list item1 item2))
                     lst2))
              lst1)))

; all subsets of a list
(define (all-subsets lst)
  (if (empty? lst)
      '(())
      (let* [(prev (all-subsets (rest lst)))]
        (append
         prev (map (lambda (item) (cons (first lst) item))
                   prev)))))

; Create a list of all assignments of #t and #f to a nonempty list of items
; (typically formulas).  Each assignment is a lookup table of entry's.

(define (all-assigns lst)
  (if (<= (length lst) 1)
      (list (list (entry (first lst) #t)) (list (entry (first lst) #f)))
      (let [(prev-assigns (all-assigns (rest lst)))]
        (append
         (map (lambda (assign) (cons (entry (first lst) #t) assign)) prev-assigns)
         (map (lambda (assign) (cons (entry (first lst) #f) assign)) prev-assigns)))))

; Sort a list of symbols according to their conversions to strings
; Example:
;> (sort-symbols '(b a aq b c))
;'(a aq b b c)

(define (sort-symbols lst)
  (sort lst
        (lambda (sym1 sym2) (string<=? (symbol->string sym1) (symbol->string sym2)))))


; is predicate true for all elements of list
(define (for-all? pred? lst)
  (cond
    [(empty? lst)
     #t]
    [(not (pred? (first lst)))
     #f]
    [else
     (for-all? pred? (rest lst))]))

; is predicate true for at least one element of list
(define (exists? pred? lst)
  (not (for-all? (lambda (item) (not (pred? item))) lst)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Various random selection methods

; pick a random element of a non-empty list
(define (pick lst)
  (list-ref lst (random (length lst))))

; Pick k distinct random elements from a list.
; List is assumed to have at least k elements.

(define (pick-k k lst)
  (if (<= k 0)
      '()
      (let ((item (pick lst)))
        (cons
         item
         (pick-k (- k 1) (remove item lst))))))

; Return a randomly chosen subset of a set

(define (random-subset lst)
  (filter (lambda (x) (<= (random) 0.5)) lst))

; Generate a sorted list of d distinct random nonnegative integers less than k.
; d is assumed to be less than k.

(define (random-sorted-list k d)
  (random-sorted-list-loop k d '()))
(define (random-sorted-list-loop k d result)
  (if (<= k 0)
      (sort result <=)
      (let [(next (random d))]
        (if (member next result)
            (random-sorted-list-loop k d result)
            (random-sorted-list-loop (- k 1) d (cons next result))))))

; Given a list of elements and their counts, select an element of the list
; with probability equal to the fraction of the total counts.

(define (random-select-counts lst counts)
  (let* [(total (apply + counts))
         (probs (map (lambda (count) (/ (* 1.0 count) total)) counts))
         (random-number (random))]
    (list-ref lst (select-index-probs random-number 0 probs))))

; Given a list of probabilities that sum to 1, select an index according to the probabilities.

(define (select-index-probs random-number index probs)
  (if (<= random-number (first probs))
      index
      (select-index-probs (- random-number (first probs)) (+ index 1) (rest probs))))
