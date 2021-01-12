#lang racket
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Definitions for strings and ultimately periodic words for automata

; Dana Angluin
; 9/27/2020

(provide all-strings-exactly-k
         all-strings-at-most-k
         bounded-seqs-k)

(require "general-utilities.rkt")

; A "string" is represented as a finite sequence of symbols.
; Example: '(a b b a b)

; Make a list of all strings of length k over given alphabet

(define (all-strings-exactly-k alphabet k)
  (if (<= k 0)
      '(())
      (let ((prev (all-strings-exactly-k alphabet (- k 1))))
        (apply append
               (map (lambda (symbol)
                      (map (lambda (str)
                             (cons symbol str))
                           prev))
                    alphabet)))))

; Make a list of all strings of length up to k over the given alphabet.

(define (all-strings-at-most-k alphabet k)
  (apply append
         (map (lambda (l) (all-strings-exactly-k alphabet l))
              (range 0 (+ 1 k)))))

; An ultimately periodic word u(v)^w is represented as a list of two "strings".
; Example: '((a b)(b b a))

; Make a list of all ultimately periodic words (u v) over the given alphabet such
; that |u| is in [0,k] and |v| is in [1,k]

(define bounded-seqs-k
  (lambda (alphabet k)
    (let* ((strings (all-strings-at-most-k alphabet k))
	   (nonempty-strings (set-diff strings '(()))))
      (c-product strings nonempty-strings))))


