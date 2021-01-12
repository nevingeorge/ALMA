#lang racket
; Procedures for SUBAs:
; Test whether an NBA is a SUBA, generate random SUBAs, write a SUBA to a file,

; September 2020
; Dana Angluin

(provide nba-is-suba? generate-random-subas write-SUBA)
 
(require "general-utilities.rkt")
(require "sequence-utilities.rkt")
(require "aut-definitions.rkt")
(require "aut-procedures.rkt")
(require "dfas-nfas-nbas.rkt")

; An example (random) SUBA with 4 states.
; It was generated using generate-random-subas with
; alphabet = '(a b c), #states = 4, #final states = 1, #transitions to remove = 4

(define ex-suba
  (aut-f
   (aut
    '(a b c)
    '(1 2 3 4)
    '(1)
    (list
     (entry '(1 a) 2)
     (entry '(1 b) 2)
     (entry '(1 c) 1)
     (entry '(2 a) 3)
     (entry '(3 a) 4)
     (entry '(3 c) 2)
     (entry '(4 b) 4)
     (entry '(4 b) 3)))
   '(3)))

; Given a dfa, determine whether it accepts the empty language.

(define (dfa-empty-lang? dfa1)
  (if (empty? (aut-f-init-states dfa1))
      #t
      (let* ((init-state (first (aut-f-init-states dfa1)))
             (reachable (aut-reachable (aut-f-aut dfa1) init-state)))
        ;    (displayln "dfa-empty-lang?")
        ;    (displayln reachable)
        (empty?
         (set-intersect reachable (aut-f-final-states dfa1))))))

; Given two dfas over the same alphabet, return #t if their languages
; have an empty intersection, and #f otherwise.  This is more
; efficient than computing the dfa for their intersection and
; testing it for emptiness.

(define (dfa-langs-empty-intersection? dfa1 dfa2)
  (if (or (empty? (aut-f-init-states dfa1))
          (empty? (aut-f-init-states dfa2)))
      #t
      (let* ((alphabet1 (aut-f-alphabet dfa1))
             (states1 (aut-f-states dfa1))
             (initial-state1 (first (aut-f-init-states dfa1)))
             (trans1 (aut-f-trans dfa1))
             (final-states1 (aut-f-final-states dfa1))
             (alphabet2 (aut-f-alphabet dfa2))
             (states2 (aut-f-states dfa2))
             (initial-state2 (first (aut-f-init-states dfa2)))
             (trans2 (aut-f-trans dfa2))
             (final-states2 (aut-f-final-states dfa2))
             (alphabet (set-intersect alphabet1 alphabet2)))
  
        (dfa-langs-empty-intersection-loop
         alphabet trans1 trans2 final-states1 final-states2
         '() (list (list initial-state1 initial-state2))))))
    
(define (dfa-langs-empty-intersection-loop
         alphabet trans1 trans2 final-states1 final-states2
         processed to-do)
  (if (or (empty? alphabet) (empty? to-do))
      #t
      (let* ((new-processed (append processed to-do))
             (next-pairs
              (remove-duplicates
               (apply append
                      (map (lambda (item)
                             (map (lambda (symbol)
                                    (let* ((result1 (lookup (list (first item)
                                                                  symbol)
                                                            trans1))
                                           (result2 (lookup (list (second item)
                                                                  symbol)
                                                            trans2)))
                                      (list result1 result2)))
                                  alphabet))
                           to-do))))
             (new-pairs
              (filter (lambda (item)
                        (and (not (equal? (first item) #f))
                             (not (equal? (second item) #f))
                             (not (member item new-processed))))
                      next-pairs))
             (final-pairs
              (filter (lambda (item)
                        (and (member (first item) final-states1)
                             (member (second item) final-states2)))
                      new-pairs)))
        
        (if (not (empty? final-pairs))
            #f
            (dfa-langs-empty-intersection-loop
             alphabet trans1 trans2 final-states1 final-states2
             new-processed new-pairs)))))

; Given two dfas over the same alphabet, return a dfa for the
; intersection of their two languages.  This shouldn't be used
; for the emptiness test, since it creates the automaton on
; all pairs of states and only then trims it.

(define (dfa-intersect dfa1 dfa2)
  (display " dfa-intersect, number of states ")
  (display (length (aut-f-states dfa1))) (display " ")
  (displayln (length (aut-f-states dfa2)))
  (if (or (empty? (aut-f-init-states dfa1))(empty? (aut-f-init-states dfa2)))
      (aut-f (aut (aut-f-alphabet dfa1) '() '() '()) '())
      (let* ((alphabet1 (aut-f-alphabet dfa1))
             (states1 (aut-f-states dfa1))
             (initial-state1 (first (aut-f-init-states dfa1)))
             (trans1 (aut-f-trans dfa1))
             (final-states1 (aut-f-final-states dfa1))
             (alphabet2 (aut-f-alphabet dfa2))
             (states2 (aut-f-states dfa2))
             (initial-state2 (first (aut-f-init-states dfa2)))
             (trans2 (aut-f-trans dfa2))
             (final-states2 (aut-f-final-states dfa2)))
        
        (let* ((alphabet (set-intersect alphabet1 alphabet2))
               (states (c-product states1 states2))
               (initial-state (list initial-state1 initial-state2))
               (final-states (c-product final-states1 final-states2))
               (trans
                (apply append
                       (map (lambda (state)
                              (map (lambda (symbol)
                                     (entry
                                      (list state symbol)
                                      (list
                                       (lookup (list (first state) symbol)
                                               trans1)
                                       (lookup (list (second state) symbol)
                                               trans2))))
                                   alphabet))
                            states))))
          (nfa-trim
           (aut-f
            (aut
             alphabet
             states
             (list initial-state)
             trans)
            final-states))))))

; Determine whether a given NBA is a SUBA.
; First decide if reverse transition is deterministic and return #f if not.
; Then check to see if there are two different states that accept u(v)^w
; for some strings u and v, with v nonempty.

(define (nba-is-suba? nba)
  (let* ((r-trans (reverse-transition (aut-f-trans nba))))
    (if (not (deterministic-trans? r-trans))
        (begin (displayln "not reverse deterministic") #f)
        (let* ((alphabet (aut-f-alphabet nba))
               (states (aut-f-states nba))
               (state-indices (range (length states)))
               (init-states (aut-f-init-states nba))
               (trans (aut-f-trans nba))
               (final-states (aut-f-final-states nba))
               (new-aut
                (aut-f
                 (aut 
                  alphabet states init-states r-trans)
                 final-states))
               (state-dfas
                (map (lambda (index)
                       (final-loop-lang new-aut (list-ref states index)))
                     state-indices))
               (all-pairs-indices (c-product state-indices state-indices))
               (distinct-pairs-indices
                (filter (lambda (item) (< (first item) (second item)))
                        all-pairs-indices)))
;          (for-each displayln state-dfas)
          (if (exists?
               (lambda (item)
;                 (display " pair: ")(display item)
                 (let* ((index1 (first item))
                        (state1 (list-ref states index1))
                        (index2 (second item))
                        (state2 (list-ref states index2))
                        (dfa1 (list-ref state-dfas index1))
                        (dfa2 (list-ref state-dfas index2)))
                   (if (not (dfa-langs-empty-intersection? dfa1 dfa2))
                       #t #f)))
               distinct-pairs-indices)
              #f
              #t)))))


; Generate a given number of "random" SUBAs, with arguments
; alphabet, number of states, number of final states, number of missing transitions.
; This generates an NBA by generating a random reverse-deterministic transition
; relation, removing a random subset of the transitions of the given size,
; and testing to see if the result is a SUBA.  The yield of actual SUBAs is
; far less than 100%.  They also typically have fewer than the given number of states.

; It prints out various messages about its progress.
; Example:
; > (define x (generate-random-subas 10 '(a b c) 30 3 30))

(define (generate-random-subas
         count alphabet number-of-states number-of-final-states number-missing)
  (generate-random-subas-loop
   count alphabet number-of-states number-of-final-states number-missing '()))

(define (generate-random-subas-loop
         count alphabet number-of-states number-of-final-states number-missing subas)
  (if (<= count 0)
      (reverse subas)
      (let* ((test-nba
              (random-reverse-det-nba-k
               alphabet number-of-states number-of-final-states number-missing)))
        (if (nba-is-suba? test-nba)
            (begin
              (generate-random-subas-loop
               (- count 1) alphabet number-of-states number-of-final-states
               number-missing (cons test-nba subas)))
            (generate-random-subas-loop
             count alphabet number-of-states number-of-final-states number-missing
             subas)))))

; Function to translate states into numbers 1 through number of states,
; as expected for input to SUBA.jar.

(define (number-of-state state list-of-states)
  (+ 1 (index-of list-of-states state)))

; Write out a SUBA to a text file in the format required by SUBA.jar.
; First line will be a comment with the given identifier.
; It also includes a comment with the initial state(s).

; Example:
; (write-SUBA "Example SUBA" ex-suba "test.txt")

(define (write-SUBA id nba1 file-name)
  (let* ((o-port (open-output-file file-name #:exists 'replace))
         (alphabet (aut-f-alphabet nba1))
         (states (aut-f-states nba1))
         (init-states (aut-f-init-states nba1))
         (trans (aut-f-trans nba1))
         (final-states (aut-f-final-states nba1)))
    (display "// id = " o-port)
    (displayln id o-port)
    (displayln "// number of states" o-port)
    (displayln (length states) o-port)
    (displayln "// alphabet" o-port)
    (for-each (lambda (symbol)
                (display symbol o-port)
                (display " " o-port))
              alphabet)
    (displayln "" o-port)
    (displayln "// initial states " o-port)
    (for-each (lambda (state)
                (display (number-of-state state states) o-port)
                (display " " o-port))
              init-states)
    (displayln "" o-port)
    (displayln "// final states" o-port)
    (for-each (lambda (state)
                (display (number-of-state state states) o-port)
                (display " " o-port))
              final-states)
    (displayln "" o-port)
    (displayln "// number of transitions" o-port)
    (displayln (length trans) o-port)
    (displayln "// transitions" o-port)
    (for-each (lambda (item)
                (let ((key (entry-key item))
                      (value (entry-value item)))
                  (display (number-of-state (first key) states) o-port)
                  (display " " o-port)
                  (display (second key) o-port)
                  (display " " o-port)
                  (displayln (number-of-state value states) o-port)))
              trans)
    (close-output-port o-port)))

(define (write-many-SUBA list-nba file-name)
  (let* ((o-port (open-output-file file-name #:exists 'replace)))
    (displayln "// number of SUBA" o-port)
    (displayln (length list-nba) o-port)
    (displayln "" o-port)
    (for-each (lambda (nba)
                (let* ((alphabet (aut-f-alphabet nba))
                       (states (aut-f-states nba))
                       (init-states (aut-f-init-states nba))
                       (trans (aut-f-trans nba))
                       (final-states (aut-f-final-states nba)))
                  (displayln "// number of states" o-port)
                  (displayln (length states) o-port)
                  (displayln "// alphabet" o-port)
                  (for-each (lambda (symbol)
                              (display symbol o-port)
                              (display " " o-port))
                            alphabet)
                  (displayln "" o-port)
                  (displayln "// final states" o-port)
                  (for-each (lambda (state)
                              (display (number-of-state state states) o-port)
                              (display " " o-port))
                            final-states)
                  (displayln "" o-port)
                  (displayln "// number of transitions" o-port)
                  (displayln (length trans) o-port)
                  (displayln "// transitions" o-port)
                  (for-each (lambda (item)
                              (let ((key (entry-key item))
                                    (value (entry-value item)))
                                (display (number-of-state (first key) states) o-port)
                                (display " " o-port)
                                (display (second key) o-port)
                                (display " " o-port)
                                (displayln (number-of-state value states) o-port)))
                            trans)
                  (displayln "" o-port))) list-nba)
    (close-output-port o-port)))

(define (write-specific-SUBA file-name)
  (let* ((o-port (open-output-file file-name #:exists 'replace))
         (list-nba '()))

    (set! list-nba (append (generate-random-subas 1000 '(a b c) 20 2 6) list-nba))
    (set! list-nba (append (generate-random-subas 1000 '(a b c) 20 2 13) list-nba))
    (set! list-nba (append (generate-random-subas 1000 '(a b c) 20 2 20) list-nba))

    (displayln "100%")

    (displayln "// number of SUBA" o-port)
    (displayln (length list-nba) o-port)
    (displayln "" o-port)
    (for-each (lambda (nba)
                (let* ((alphabet (aut-f-alphabet nba))
                       (states (aut-f-states nba))
                       (init-states (aut-f-init-states nba))
                       (trans (aut-f-trans nba))
                       (final-states (aut-f-final-states nba)))
                  (displayln "// number of states" o-port)
                  (displayln (length states) o-port)
                  (displayln "// alphabet" o-port)
                  (for-each (lambda (symbol)
                              (display symbol o-port)
                              (display " " o-port))
                            alphabet)
                  (displayln "" o-port)
                  (displayln "// final states" o-port)
                  (for-each (lambda (state)
                              (display (number-of-state state states) o-port)
                              (display " " o-port))
                            final-states)
                  (displayln "" o-port)
                  (displayln "// number of transitions" o-port)
                  (displayln (length trans) o-port)
                  (displayln "// transitions" o-port)
                  (for-each (lambda (item)
                              (let ((key (entry-key item))
                                    (value (entry-value item)))
                                (display (number-of-state (first key) states) o-port)
                                (display " " o-port)
                                (display (second key) o-port)
                                (display " " o-port)
                                (displayln (number-of-state value states) o-port)))
                            trans)
                  (displayln "" o-port))) list-nba)
    (close-output-port o-port)))