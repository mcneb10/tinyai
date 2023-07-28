pipeline {
    agent { label 'linux' }
    stages {
        stage('Setup') {
            steps {
                echo 'Setting up...'
                echo '- Updating'
                sh 'sudo apt-get update'
                echo '- Installing build tools'
                sh 'sudo apt-get install -y build-essential'
            }
        }
        stage('Build') {
            steps {
                echo 'Building...'
                sh 'make -C experiments/xor/'
            }
        }
        stage('Test') {
            steps {
                echo 'Testing...'
                sh './experiments/xor/xor'
                sh './experiments/xor/xor_independent'
                sh './experiments/xor/xor_changeable_activation_aggregation'
                sh './experiments/xor/xor_changeable_activation_aggregation_independent'
            }
        }
    }
}